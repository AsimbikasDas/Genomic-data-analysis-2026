import pandas as pd
import wandb
import torch
import os
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split

from models.resnet_1d import RNAResNet1D
from models.bilstm import RNABiLSTM
from models.oned_cnn import RNACnn

def get_device():
    if torch.cuda.is_available():
        return torch.device("cuda")
    elif torch.backends.mps.is_available():
        return torch.device("mps")
    else:
        return torch.device("cpu")

class GenomeDataset(Dataset):
    def __init__(self, sequences, labels):

        self.char_to_int = {
            'A': 0,
            'C': 1,
            'G': 2,
            'T': 3,
            'U': 3,
            'N': 4,
        }

        print("Converting sequences to integers...")
        self.encoded_sequences = [
            [self.char_to_int.get(char, 4) for char in seq]
            for seq in sequences
        ]

        self.seq_tensors = torch.tensor(self.encoded_sequences, dtype=torch.long)
        self.labels = torch.tensor(labels, dtype=torch.float32)

    def __len__(self):
        return len(self.labels)
    
    def __getitem__(self, idx):
        return self.seq_tensors[idx], self.labels[idx].unsqueeze(0)

if __name__=="__main__":

    BATCH_SIZE = 64
    learning_rate = 0.001
    EPOCHS = 50
    threshold = 0.45
    patience = 7
    imbalance_penalty = 4.0
    architecture = "CNN" # Options: "CNN", "ResNet", "BiLSTM"

    device = get_device()

    wandb.init(
        project=f"rna-{architecture}",
        config={
            "epochs": EPOCHS,
            "batch_size": BATCH_SIZE,
            "learning_rate": learning_rate,
            "architecture": architecture,
            "imbalance_penalty": imbalance_penalty
        }
    )

    int("Loading data...")
    data_path = "data/combined_exon_mirna_2.csv"
    df = pd.read_csv(data_path)

    df['label'] = df['label'].map({'exon': 0, 'miRNA': 1})

    X_train, X_test, y_train, y_test = train_test_split(
        df['sequence'].values,
        df['label'].values,
        test_size=0.2,
        random_state=42,
        stratify=df['label'].values
    )

    train_dataset = GenomeDataset(X_train, y_train)
    test_dataset = GenomeDataset(X_test, y_test)

    train_dataloader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True, num_workers=4)
    test_dataloader = DataLoader(test_dataset, batch_size=BATCH_SIZE, shuffle=False, num_workers=4)
    
    if architecture == "CNN":
        model = RNACnn().to(device)
    elif architecture == "ResNet":
        model = RNAResNet1D().to(device)
    else:
        model = RNABiLSTM().to(device)

    pos_weight = torch.tensor([imbalance_penalty], device=device)
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    optimizer = optim.AdamW(model.parameters(), lr=learning_rate, weight_decay=1e-4)

    best_val_loss = float('inf')
    best_recall_score = 0.0
    epochs_no_improve = 0

    print("Starting training...")
    for epoch in range(EPOCHS):
        model.train()
        running_loss = 0

        for i, (inputs, labels) in enumerate(train_dataloader):
            inputs, labels = inputs.to(device), labels.to(device)

            inputs = F.one_hot(inputs, num_classes=5).float()
            inputs = inputs[:, :, :4]
            if architecture != "BiLSTM":
                inputs = inputs.transpose(1, 2).contiguous()

            optimizer.zero_grad()

            outputs = model(inputs)

            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()

            running_loss += loss.item()

        avg_loss = running_loss / len(train_dataloader)
        print(f"Epoch [{epoch+1}/{EPOCHS}], Training Loss: {avg_loss:.4f}")

        wandb.log({
            "training/loss": avg_loss,
            "epoch": epoch
        })

        model.eval()
        
        val_loss = 0.0
        tp, tn, fp, fn = 0, 0, 0, 0

        all_probs = []
        all_labels = []

        with torch.no_grad():

            for inputs, labels in test_dataloader:
                inputs, labels = inputs.to(device), labels.to(device)

                inputs = F.one_hot(inputs, num_classes=5).float()
                inputs = inputs[:, :, :4]
                if architecture != "BiLSTM":
                    inputs = inputs.transpose(1, 2).contiguous()

                outputs = model(inputs)
                loss = criterion(outputs, labels)
                val_loss += loss.item()

                probs = torch.sigmoid(outputs)
                all_probs.extend(probs.cpu().numpy().flatten())
                all_labels.extend(labels.cpu().numpy().flatten())

                predicted = (probs > threshold).float()

                tp += ((predicted == 1) & (labels == 1)).sum().item()
                tn += ((predicted == 0) & (labels == 0)).sum().item()
                fp += ((predicted == 1) & (labels == 0)).sum().item()
                fn += ((predicted == 0) & (labels == 1)).sum().item()

        avg_loss = val_loss / len(test_dataloader)
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        print(f"Epoch [{epoch+1}/{EPOCHS}], Validation Precision: {precision:.2f}, Recall: {recall:.2f}, F1: {f1:.2f}")


        wandb.log({
            "validation/precision": precision,
            "validation/recall": recall,
            "validation/f1": f1,
            "pr_curve": wandb.plot.pr_curve(y_true=all_labels, y_probas=[[1-p, p] for p in all_probs], labels=["Exon", "miRNA"]),
            "epoch": epoch
        })

        if recall > best_recall_score:
            best_recall_score = recall
            print(f"New best validation recall score: {best_recall_score:.4f}.")
            model_path = "models/cnn_1d/recall_max/best_model_cnn_5.pth"
            os.makedirs(os.path.dirname(model_path), exist_ok=True)
            torch.save(model.state_dict(), model_path)
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1
            print(f"No improvement in validation recall score for {epochs_no_improve} epoch(s).")
            if epochs_no_improve >= patience:
                print("Early stopping triggered.")
                break

    print("Training complete!")        
