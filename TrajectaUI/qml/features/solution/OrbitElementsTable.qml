pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: table
    required property SolutionViewModel solution
    readonly property OrbitModel orbit: solution.orbit
    property bool equinoctial: false
    readonly property bool equinoctialAvailable: solution.equinoctialElements.length === 6
    onEquinoctialAvailableChanged: if (!equinoctialAvailable) equinoctial = false
    spacing: 0
    SectionHeading { title: qsTr("Nominal orbit"); Layout.fillWidth: true; Layout.bottomMargin: 14 }
    CustomSegmentedControl {
        objectName: "orbit-elements-mode"
        options: [table.solution.cometary ? qsTr("Cometary") : qsTr("Keplerian"), qsTr("Equinoctial")]
        disabledIndices: table.equinoctialAvailable ? [] : [1]
        currentIndex: table.equinoctial ? 1 : 0
        Layout.bottomMargin: 9
        onActivated: function(index) { table.equinoctial = index === 1 && table.equinoctialAvailable }
    }
    RowLayout {
        Layout.fillWidth: true
        Layout.preferredHeight: 34
        spacing: 10
        CustomText { text: qsTr("Parameter"); secondary: true; font.pixelSize: 10; Layout.preferredWidth: 50 }
        CustomText { text: qsTr("Value"); secondary: true; font.pixelSize: 10; Layout.fillWidth: true; horizontalAlignment: Text.AlignRight }
        CustomText { text: "σ"; secondary: true; font.pixelSize: 10; Layout.preferredWidth: 72; horizontalAlignment: Text.AlignRight }
    }
    Repeater {
        objectName: "orbit-element-rows"
        model: table.equinoctial ? table.solution.equinoctialElements : [
            {label: table.solution.cometary ? "q" : "a", unit: "AU", value: table.solution.cometary ? table.orbit.perigee : table.orbit.semiMajorAxis, sigma: table.solution.cometary ? table.orbit.sigmaQ : table.orbit.sigmaA},
            {label: "e", unit: "", value: table.orbit.eccentricity, sigma: table.orbit.sigmaE},
            {label: "i", unit: "°", value: table.orbit.inclination, sigma: table.orbit.sigmaI},
            {label: "Ω", unit: "°", value: table.orbit.ascendingNode, sigma: table.orbit.sigmaNode},
            {label: "ω", unit: "°", value: table.orbit.argPericentre, sigma: table.orbit.sigmaW},
            {label: table.solution.cometary ? "ν" : "M", unit: "°", value: table.solution.cometary ? table.orbit.trueAnomaly : table.orbit.meanAnomaly, sigma: table.solution.cometary ? table.orbit.sigmaNu : table.orbit.sigmaM}
        ]
        Item {
            id: elementRow
            required property var modelData
            required property int index
            objectName: "orbit-element-" + index
            Layout.fillWidth: true
            implicitHeight: 42
            Rectangle { width: parent.width; height: 1; color: AppTheme.lineSoft }
            RowLayout {
                anchors.fill: parent
                spacing: 10
                RowLayout {
                    Layout.preferredWidth: 50
                    spacing: 4
                    CustomText { text: elementRow.modelData.label; font.family: Qt.platform.os === "windows" ? "Georgia" : UI.font; font.italic: true; font.pixelSize: 14 }
                    CustomText { text: elementRow.modelData.unit; visible: text.length > 0; secondary: true; font.pixelSize: 10 }
                    Item { Layout.fillWidth: true }
                }
                CustomText { text: UI.number(elementRow.modelData.value, 9); numeric: true; font.pixelSize: 11; Layout.fillWidth: true; horizontalAlignment: Text.AlignRight }
                CustomText { text: table.equinoctial || table.solution.hasSigmas ? UI.sigma(elementRow.modelData.sigma) : "—"; numeric: true; secondary: true; font.pixelSize: 11; Layout.preferredWidth: 72; horizontalAlignment: Text.AlignRight }
            }
        }
    }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.lineSoft }
    Item {
        Layout.fillWidth: true
        implicitHeight: 42
        visible: table.orbit.eccentricity < 1 && Number.isFinite(table.orbit.orbitalPeriod) && table.orbit.orbitalPeriod > 0
        RowLayout {
            anchors.fill: parent
            spacing: 10
            CustomText { text: qsTr("Period"); font.pixelSize: 11; Layout.preferredWidth: 50 }
            CustomText {
                text: table.orbit.orbitalPeriod > 365.25
                      ? qsTr("%1 days (%2 years)").arg(UI.number(table.orbit.orbitalPeriod, 4)).arg(UI.number(table.orbit.orbitalPeriod / 365.25, 4))
                      : qsTr("%1 days").arg(UI.number(table.orbit.orbitalPeriod, 4))
                numeric: true
                font.pixelSize: 11
                Layout.fillWidth: true
                horizontalAlignment: Text.AlignRight
            }
            CustomText { text: table.solution.hasSigmas ? UI.sigma(table.orbit.sigmaPeriod) : "—"; numeric: true; secondary: true; font.pixelSize: 11; Layout.preferredWidth: 72; horizontalAlignment: Text.AlignRight }
        }
    }
    Rectangle {
        Layout.fillWidth: true
        implicitHeight: 1
        color: AppTheme.lineSoft
        visible: table.orbit.eccentricity < 1 && Number.isFinite(table.orbit.orbitalPeriod) && table.orbit.orbitalPeriod > 0
    }
    NonGravParameters {
        objectName: "nominal-non-grav"
        solution: table.solution
        Layout.fillWidth: true
        Layout.topMargin: 16
        Layout.bottomMargin: 16
    }
    CustomTextLink { text: qsTr("All parameters and vectors"); implicitHeight: 16; Layout.topMargin: 12; onClicked: details.open() }
    MoidSummary { solution: table.solution; Layout.fillWidth: true; Layout.topMargin: 17 }
    OrbitDetailsDialog { id: details; solution: table.solution }
}
