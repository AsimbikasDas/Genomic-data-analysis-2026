import torch.nn as nn
import torch.nn.functional as F

class RNACnn(nn.Module):
    def __init__(self):
        super(RNACnn, self).__init__()

        self.conv1 = nn.Conv1d(in_channels=4, out_channels=128, kernel_size=4)

        self.pool = nn.MaxPool1d(kernel_size=2)

        self.conv2 = nn.Conv1d(in_channels=128, out_channels=256, kernel_size=3)

        self.adaptive_pool = nn.AdaptiveAvgPool1d(1)

        self.flatten = nn.Flatten()

        self.fc1 = nn.Linear(256, 64)
        self.dropout = nn.Dropout(0.3)
        self.fc2 = nn.Linear(64, 1)

    def forward(self, x):
        x = F.relu(self.conv1(x))
        x = self.pool(x)
        x = F.relu(self.conv2(x))
        x = self.adaptive_pool(x)
        x = self.flatten(x)
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.fc2(x)
        return x
