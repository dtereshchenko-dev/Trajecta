pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

Flickable {
    id: panel
    objectName: "covariance-panel"
    required property CovarianceTableModel matrix
    readonly property bool stacked: width < 1040
    contentWidth: width
    contentHeight: content.implicitHeight + 24
    clip: true
    boundsBehavior: Flickable.StopAtBounds
    ScrollBar.vertical: CustomScrollBar {}

    ColumnLayout {
        id: content
        width: panel.width
        spacing: 0
        RowLayout {
            Layout.fillWidth: true
            spacing: 16
            CustomText { text: qsTr("Solution covariance"); font.pixelSize: 17; font.weight: Font.DemiBold; Layout.fillWidth: true }
            CustomSegmentedControl {
                objectName: "covariance-mode"
                options: [qsTr("Correlation"), qsTr("COV")]
                currentIndex: panel.matrix.correlation ? 0 : 1
                enabled: panel.matrix.dimension > 0
                onActivated: function(index) { panel.matrix.correlation = index === 0 }
            }
        }
        Flow {
            Layout.fillWidth: true
            Layout.topMargin: 12
            Layout.bottomMargin: 22
            spacing: 14
            CustomComboBox {
                objectName: "covariance-basis"
                width: 190
                model: panel.matrix.equinoctialAvailable ? [qsTr("Equinoctial elements"), qsTr("Cartesian state")] : [qsTr("Cartesian state")]
                currentIndex: panel.matrix.equinoctialAvailable && !panel.matrix.equinoctial ? 1 : 0
                enabled: panel.matrix.equinoctialAvailable
                onActivated: function(index) { panel.matrix.equinoctial = index === 0 && panel.matrix.equinoctialAvailable }
            }
            InlineValue { height: 33; label: qsTr("Epoch"); value: panel.matrix.epoch || "\u2014" }
            InlineValue { height: 33; label: qsTr("COV scale"); value: Number.isFinite(panel.matrix.objectiveScale) ? String(panel.matrix.objectiveScale) : "\u2014" }
        }
        CustomText {
            visible: panel.matrix.dimension === 0
            text: qsTr("No covariance available")
            secondary: true
            Layout.fillWidth: true
            Layout.topMargin: 32
        }
        GridLayout {
            Layout.fillWidth: true
            visible: panel.matrix.dimension > 0
            columns: panel.stacked ? 1 : 2
            columnSpacing: 36
            rowSpacing: 24
            ColumnLayout {
                Layout.fillWidth: true
                Layout.preferredWidth: panel.stacked ? panel.width : (panel.width - 36) * 1.7 / 2.7
                Layout.minimumWidth: 0
                Layout.alignment: Qt.AlignTop
                spacing: 0
                CovarianceMatrix { Layout.fillWidth: true; matrix: panel.matrix }
                RowLayout {
                    Layout.fillWidth: true
                    Layout.topMargin: 12
                    PlotLegendItem { label: "\u22121.0"; color: AppTheme.blue }
                    Item { Layout.fillWidth: true }
                    CustomText { text: "0"; secondary: true; font.pixelSize: 11 }
                    Item { Layout.fillWidth: true }
                    PlotLegendItem { label: "+1.0"; color: AppTheme.accent }
                }
                Item {
                    Layout.fillWidth: true
                    Layout.topMargin: 28
                    implicitHeight: units.implicitHeight + 24
                    Rectangle { width: 2; height: parent.height; color: AppTheme.muted }
                    CustomText {
                        id: units
                        objectName: "covariance-units"
                        x: 16
                        y: 12
                        width: parent.width - 30
                        text: panel.matrix.basisUnits
                        secondary: true
                        font.pixelSize: 11
                        wrapMode: Text.Wrap
                        elide: Text.ElideNone
                    }
                }
            }
            Item {
                objectName: "covariance-inspector"
                Layout.fillWidth: true
                Layout.preferredWidth: panel.stacked ? panel.width : (panel.width - 36) / 2.7
                Layout.minimumWidth: panel.stacked ? 0 : 260
                Layout.fillHeight: !panel.stacked
                Layout.alignment: Qt.AlignTop
                implicitHeight: inspector.implicitHeight + (panel.stacked ? 24 : 0)
                Rectangle { width: panel.stacked ? parent.width : 1; height: panel.stacked ? 1 : parent.height; color: AppTheme.line }
                CovarianceDetails {
                    id: inspector
                    x: panel.stacked ? 0 : 24
                    y: panel.stacked ? 24 : 0
                    width: parent.width - x
                    matrix: panel.matrix
                }
            }
        }
    }
}
