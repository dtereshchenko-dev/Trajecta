import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

BaseDialog {
    id: dialog
    objectName: "orbit-details-dialog"
    required property SolutionViewModel solution
    title: qsTr("Nominal orbital parameters")
    acceptText: qsTr("Close")
    width: Math.min(1000, parent ? parent.width - 56 : 1000)
    height: Math.min(Math.ceil(implicitHeight), parent ? parent.height - 56 : Math.ceil(implicitHeight))
    contentItem: ScrollView {
        id: scroll
        contentWidth: availableWidth
        implicitHeight: content.implicitHeight
        rightPadding: 14
        clip: true
        ScrollBar.vertical: CustomScrollBar {
            objectName: "orbit-details-scrollbar"
            parent: scroll
            x: scroll.width - width
            y: scroll.topPadding
            height: scroll.availableHeight
            policy: scroll.contentHeight > scroll.availableHeight + 1 ? ScrollBar.AlwaysOn : ScrollBar.AsNeeded
        }
        ColumnLayout {
            id: content
            width: scroll.availableWidth
            spacing: 24
            RowLayout {
                Layout.fillWidth: true
                spacing: 16
                InlineValue { label: qsTr("Epoch"); value: dialog.solution.epochMjd; unit: "MJD TDB" }
                CustomText { text: "ECLM J2000"; secondary: true; font.pixelSize: 10; Layout.fillWidth: true; horizontalAlignment: Text.AlignRight }
            }
            GridLayout {
                Layout.fillWidth: true
                columns: width < 720 ? 1 : 2
                columnSpacing: 28
                rowSpacing: 24
                uniformCellWidths: true
                DerivedOrbitParameters { orbit: dialog.solution.orbit; Layout.fillWidth: true; Layout.alignment: Qt.AlignTop }
                CartesianStateTable { values: dialog.solution.stateVector; Layout.fillWidth: true; Layout.alignment: Qt.AlignTop }
            }
            ColumnLayout {
                Layout.fillWidth: true
                spacing: 12
                SectionHeading { title: qsTr("Photometry"); Layout.fillWidth: true }
                ParameterTable {
                    objectName: "photometry-parameters"
                    Layout.fillWidth: true
                    rowPrefix: "photometry"
                    unitShare: Math.max(0.08, 56 / Math.max(1, columnSpace))
                    rows: dialog.solution.cometPhotometry ? [
                        {label: "M1", value: UI.number(dialog.solution.cometTotalMagnitude, 3), unit: "mag"},
                        {label: "M2", value: UI.number(dialog.solution.cometNuclearMagnitude, 3), unit: "mag"},
                        {label: "K1", value: UI.number(dialog.solution.cometK1, 2), unit: qsTr("fixed")},
                        {label: "K2", value: UI.number(dialog.solution.cometK2, 2), unit: qsTr("fixed")},
                        {label: "PHCOF", value: UI.number(dialog.solution.cometNuclearPhaseCoefficient, 3), unit: "mag/deg"}
                    ] : [
                        {label: "H", value: UI.number(dialog.solution.absoluteMagnitude, 3), unit: "mag"},
                        {label: "G", value: UI.number(dialog.solution.asteroidSlopeG, 3), unit: qsTr("fixed")}
                    ]
                }
            }
            Rectangle { visible: estimatedSize.visible; Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
            EstimatedSize { id: estimatedSize; solution: dialog.solution; Layout.fillWidth: true }
            Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
            NonGravParameters {
                objectName: "detailed-non-grav"
                solution: dialog.solution
                detailed: true
                Layout.fillWidth: true
            }
        }
    }
}
