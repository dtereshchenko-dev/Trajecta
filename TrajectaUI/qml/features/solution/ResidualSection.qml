import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: section
    required property SolutionViewModel solution
    property real viewportWidth: 1440
    property bool timeMode: false
    spacing: 0
    RowLayout {
        Layout.fillWidth: true
        Layout.bottomMargin: 14
        spacing: 12
        CustomText { text: qsTr("Observation residuals"); font.pixelSize: 14; font.weight: Font.DemiBold; Layout.fillWidth: true }
        CustomSegmentedControl { objectName: "residual-mode"; options: [qsTr("Scatter"), qsTr("By time")]; currentIndex: section.timeMode ? 1 : 0; onActivated: function(index) { section.timeMode = index === 1 } }
    }
    Rectangle {
        Layout.fillWidth: true
        implicitHeight: plot.height + 34
        radius: 4
        color: AppTheme.surface
        border.color: AppTheme.line
        RowLayout {
            x: 13; y: 1; height: 32; width: parent.width - 26
            CustomText { text: section.timeMode ? "Δα cos δ, Δδ / UTC" : "Δα cos δ / Δδ"; secondary: true; font.pixelSize: 10 }
            Item { Layout.fillWidth: true }
            CustomText { visible: !section.timeMode && section.solution.rms > 0; text: UI.number(section.solution.rms, 2) + "″ RMS"; color: AppTheme.blue; font.pixelSize: 10 }
            CustomText { visible: !section.timeMode && section.solution.rms > 0; text: UI.number(section.solution.rms * 3, 2) + "″ 3×RMS"; color: AppTheme.amber; font.pixelSize: 10 }
            CustomText { text: qsTr("arcsec"); secondary: true; font.pixelSize: 10 }
        }
        Rectangle { x: 1; y: 32; height: 1; width: parent.width - 2; color: AppTheme.lineSoft }
        ResidualPlot {
            id: plot
            x: 1; y: 33; width: parent.width - 2
            height: section.viewportWidth >= 1700 ? 420 : section.viewportWidth <= 1200 ? 325 : 358
            observations: section.solution.observations
            referenceRms: section.solution.rms
            timeMode: section.timeMode
        }
    }
    RowLayout {
        Layout.fillWidth: true
        Layout.topMargin: 11
        Layout.bottomMargin: 11
        spacing: 16
        PlotLegendItem { label: section.timeMode ? "Δα cos δ" : qsTr("Included"); color: AppTheme.accent }
        PlotLegendItem { visible: section.timeMode; label: "Δδ"; color: AppTheme.blue }
        PlotLegendItem { label: qsTr("Excluded"); color: AppTheme.red }
        Item { Layout.fillWidth: true }
        CustomText { text: qsTr("Selected: %1").arg(section.solution.observations.selectedSourceRow < 0 ? "—" : section.solution.observations.selectedSourceRow + 1); secondary: true; font.pixelSize: 11 }
    }
}
