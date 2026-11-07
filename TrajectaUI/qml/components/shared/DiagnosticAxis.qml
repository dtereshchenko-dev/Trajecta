import QtQuick
import QtCharts
import TrajectaUI
import "../../shared/UI.js" as UI

ValueAxis {
    tickCount: 5
    labelsColor: AppTheme.chartTick
    titleBrush: AppTheme.muted
    gridLineColor: AppTheme.chartGrid
    color: AppTheme.chartGrid
    labelsFont.family: Qt.platform.os === "windows" ? "Consolas" : UI.monoFont
    labelsFont.pixelSize: 10
    titleFont.family: Qt.platform.os === "windows" ? "Segoe UI" : UI.font
    titleFont.pixelSize: 10
    labelFormat: "%.2f"
}
