import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

BaseDialog {
    id: dialog
    objectName: "orbit-details-dialog"
    required property SolutionViewModel solution
    required property var elementRows
    property bool equinoctial: false
    title: qsTr("Nominal orbital parameters")
    acceptText: qsTr("Close")
    width: Math.min(1000, parent ? parent.width - 56 : 1000)
    height: Math.min(Math.ceil(implicitHeight), parent ? parent.height - 56 : Math.ceil(implicitHeight))
    contentItem: Flickable {
        id: scroll
        contentWidth: width
        contentHeight: content.implicitHeight
        implicitHeight: content.implicitHeight
        interactive: true
        flickableDirection: Flickable.VerticalFlick
        boundsBehavior: Flickable.StopAtBounds
        clip: true
        ScrollBar.vertical: CustomScrollBar {
            objectName: "orbit-details-scrollbar"
            interactive: true
            policy: scroll.contentHeight > scroll.height ? ScrollBar.AlwaysOn : ScrollBar.AlwaysOff
        }
        ColumnLayout {
            id: content
            width: Math.max(0, scroll.width - 14)
            spacing: 24
            RowLayout {
                Layout.fillWidth: true
                spacing: 16
                InlineValue { label: qsTr("Epoch"); value: dialog.solution.epochMjd; unit: "MJD TDB" }
                CustomText { text: "ECLM J2000"; secondary: true; font.pixelSize: 10; Layout.fillWidth: true; horizontalAlignment: Text.AlignRight }
            }
            ColumnLayout {
                Layout.fillWidth: true
                spacing: 12
                SectionHeading { title: dialog.equinoctial ? qsTr("Equinoctial elements") : dialog.solution.cometary ? qsTr("Cometary elements") : qsTr("Keplerian elements"); Layout.fillWidth: true }
                ParameterTable {
                    objectName: "full-precision-elements"
                    Layout.fillWidth: true
                    rowPrefix: "full-element"
                    showSigma: true
                    selectableValues: true
                    rows: dialog.elementRows.map(row => ({
                        label: row.label, value: UI.fullPrecision(row.value), unit: row.unit,
                        sigma: dialog.equinoctial || dialog.solution.hasSigmas ? UI.sigma(row.sigma) : "\u2014"
                    }))
                }
            }
            DerivedOrbitParameters { orbit: dialog.solution.orbit; hasSigmas: dialog.solution.hasSigmas; Layout.fillWidth: true }
            CartesianStateTable { values: dialog.solution.stateVector; sigmas: dialog.solution.stateVectorSigmas; Layout.fillWidth: true }
            ColumnLayout {
                Layout.fillWidth: true
                spacing: 12
                SectionHeading { title: qsTr("Photometry"); Layout.fillWidth: true }
                ParameterTable {
                    objectName: "photometry-parameters"
                    Layout.fillWidth: true
                    rowPrefix: "photometry"
                    showSigma: true
                    showOrbitSigma: true
                    rows: dialog.solution.cometPhotometry ? [
                        {label: "M1", value: UI.number(dialog.solution.cometTotalMagnitude, 3), sigma: UI.sigma(dialog.solution.cometTotalMagnitudePhotometricSigma), orbitSigma: UI.sigma(dialog.solution.cometTotalMagnitudeSigma), unit: "mag"},
                        {label: "M2", value: UI.number(dialog.solution.cometNuclearMagnitude, 3), sigma: UI.sigma(dialog.solution.cometNuclearMagnitudePhotometricSigma), orbitSigma: UI.sigma(dialog.solution.cometNuclearMagnitudeSigma), unit: "mag"},
                        {label: "K1", value: UI.number(dialog.solution.cometK1, 2), sigma: "\u2014", orbitSigma: "\u2014", unit: qsTr("fixed")},
                        {label: "K2", value: UI.number(dialog.solution.cometK2, 2), sigma: "\u2014", orbitSigma: "\u2014", unit: qsTr("fixed")},
                        {label: "PHCOF", value: UI.number(dialog.solution.cometNuclearPhaseCoefficient, 3), sigma: "\u2014", orbitSigma: "\u2014", unit: "mag/deg"}
                    ] : [
                        {label: "H", value: UI.number(dialog.solution.absoluteMagnitude, 3), sigma: UI.sigma(dialog.solution.absoluteMagnitudePhotometricSigma), orbitSigma: UI.sigma(dialog.solution.absoluteMagnitudeSigma), unit: "mag"},
                        {label: "G", value: UI.number(dialog.solution.asteroidSlopeG, 3), sigma: "\u2014", orbitSigma: "\u2014", unit: qsTr("fixed")}
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
