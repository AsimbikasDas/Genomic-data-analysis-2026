import torch.nn as nn

class ResidualBlock1D(nn.Module):
    def __init__(self, in_channels, out_channels, stride=1, dropout_rate=0.2):
        super(ResidualBlock1D, self).__init__()

        self.conv1 = nn.Conv1d(
            in_channels=in_channels,
            out_channels=out_channels,
            kernel_size=3,
            stride=stride,
            padding=1,
            bias=False
        )
        self.bn1 = nn.BatchNorm1d(out_channels)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(dropout_rate)

        self.conv2 = nn.Conv1d(
            in_channels=out_channels,
            out_channels=out_channels,
            kernel_size=3,
            stride=1,
            padding=1,
            bias=False
        )
        self.bn2 = nn.BatchNorm1d(out_channels)

        self.shortcut = nn.Sequential()
        if stride != 1 or in_channels != out_channels:
            self.shortcut = nn.Sequential(
                nn.Conv1d(in_channels=in_channels, out_channels=out_channels, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm1d(out_channels)
            )

    def forward(self, x):
        identity = self.shortcut(x)

        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)
        out = self.dropout(out)

        out = self.conv2(out)
        out = self.bn2(out)
        
        out += identity
        out = self.relu(out)

        return out
    
class RNAResNet1D(nn.Module):
    def __init__(self):
        super(RNAResNet1D, self).__init__()

        self.conv1 = nn.Conv1d(
            in_channels=4,
            out_channels=64,
            kernel_size=5,
            stride=1,
            padding=2,
            bias=False
        )
        self.bn1 = nn.BatchNorm1d(64)
        self.relu = nn.ReLU()

        self.layer1 = ResidualBlock1D(in_channels=64, out_channels=64, stride=1)
        self.layer2 = ResidualBlock1D(in_channels=64, out_channels=128, stride=2)
        self.layer3 = ResidualBlock1D(in_channels=128, out_channels=256, stride=2)

        self.global_pool = nn.AdaptiveAvgPool1d(1)
        self.fc = nn.Linear(256, 1)

    def forward(self, x):
        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu(x)

        x = self.layer1(x)
        x = self.layer2(x)
        x = self.layer3(x)

        x = self.global_pool(x)
        x = x.view(x.size(0), -1)
        x = self.fc(x)

        return x
