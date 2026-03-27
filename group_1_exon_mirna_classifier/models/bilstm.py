import torch.nn as nn
import torch
import torch.nn.functional as F

class RNABiLSTM(nn.Module):
    def __init__(self, input_size=4, hidden_size=128, num_layers=2, dropout=0.3):
        super(RNABiLSTM, self).__init__()

        self.hidden_size = hidden_size
        self.num_layers = num_layers

        self.lstm = nn.LSTM(
            input_size=input_size,
            hidden_size=hidden_size,
            num_layers=num_layers,
            batch_first=True,
            bidirectional=True,
            dropout=dropout if num_layers > 1 else 0
        )

        self.fc1 = nn.Linear(hidden_size*2, 32)
        self.dropout = nn.Dropout(dropout)
        self.fc2 = nn.Linear(32, 1)

    def forward(self, x):
        _, (h_n, _) = self.lstm(x)

        hidden_forward = h_n[-2,:,:]
        hidden_backward = h_n[-1,:,:]

        final_hidden = torch.cat((hidden_forward, hidden_backward), dim=1)

        x = F.relu(self.fc1(final_hidden))
        x = self.dropout(x)
        x = self.fc2(x)
        return x
