import QtQuick
import QtCharts
import TrajectaUI

ValueAxis {
    tickCount: 5
    labelsColor: AppTheme.chartTick
    titleBrush: AppTheme.muted
    gridLineColor: AppTheme.chartGrid
    color: AppTheme.chartGrid
    labelsFont.family: "Consolas"
    labelsFont.pixelSize: 10
    titleFont.family: "Segoe UI"
    titleFont.pixelSize: 10
    labelFormat: "%.2f"
}
