pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: table
    required property SolutionViewModel solution
    readonly property OrbitModel orbit: solution.orbit
    property bool equinoctial: false
    readonly property bool equinoctialAvailable: solution.equinoctialElements.length === 6
    readonly property var elementRows: equinoctial ? solution.equinoctialElements : [
        {label: solution.cometary ? "q" : "a", unit: "AU", value: solution.cometary ? orbit.perigee : orbit.semiMajorAxis, sigma: solution.cometary ? orbit.sigmaQ : orbit.sigmaA},
        {label: "e", unit: "", value: orbit.eccentricity, sigma: orbit.sigmaE},
        {label: "i", unit: "°", value: orbit.inclination, sigma: orbit.sigmaI},
        {label: "Ω", unit: "°", value: orbit.ascendingNode, sigma: orbit.sigmaNode},
        {label: "ω", unit: "°", value: orbit.argPericentre, sigma: orbit.sigmaW},
        {label: solution.cometary ? "ν" : "M", unit: "°", value: solution.cometary ? orbit.trueAnomaly : orbit.meanAnomaly, sigma: solution.cometary ? orbit.sigmaNu : orbit.sigmaM}
    ]
    onEquinoctialAvailableChanged: if (!equinoctialAvailable) equinoctial = false
    spacing: 0
    SectionHeading {
        title: qsTr("Nominal orbit")
        detail: qsTr("Epoch: %1 MJD TDB").arg(UI.number(parseFloat(table.solution.epochMjd), 1))
        detailFontSize: 12
        Layout.fillWidth: true
        Layout.bottomMargin: 12
    }
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
        CustomText { text: qsTr("Parameter"); secondary: true; font.pixelSize: 10; Layout.preferredWidth: 92 }
        CustomText { text: qsTr("Value"); secondary: true; font.pixelSize: 10; Layout.fillWidth: true; horizontalAlignment: Text.AlignRight }
        CustomText { text: "σ"; secondary: true; font.pixelSize: 10; Layout.preferredWidth: 96; horizontalAlignment: Text.AlignRight }
    }
    Repeater {
        objectName: "orbit-element-rows"
        model: table.elementRows
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
                    Layout.preferredWidth: 92
                    spacing: 4
                    CustomText { text: elementRow.modelData.label; font.family: Qt.platform.os === "windows" ? "Georgia" : UI.font; font.italic: true; font.pixelSize: 14 }
                    CustomText { text: elementRow.modelData.unit; visible: text.length > 0; secondary: true; font.pixelSize: 10 }
                    Item { Layout.fillWidth: true }
                }
                CustomText { text: UI.number(elementRow.modelData.value, 9); numeric: true; font.pixelSize: 11; Layout.fillWidth: true; horizontalAlignment: Text.AlignRight }
                CustomText { text: table.equinoctial || table.solution.hasSigmas ? UI.sigma(elementRow.modelData.sigma) : "—"; numeric: true; secondary: true; font.pixelSize: 11; Layout.preferredWidth: 96; horizontalAlignment: Text.AlignRight }
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
            CustomText { text: qsTr("Period"); font.pixelSize: 11; Layout.preferredWidth: 92 }
            CustomText {
                text: qsTr("%1 days (%2 years)").arg(UI.number(table.orbit.orbitalPeriod, 4)).arg(UI.number(table.orbit.orbitalPeriod / 365.25, 4))
                numeric: true
                font.pixelSize: 11
                Layout.fillWidth: true
                Layout.minimumWidth: 0
                horizontalAlignment: Text.AlignRight
                wrapMode: Text.Wrap
                elide: Text.ElideNone
            }
            CustomText {
                text: table.solution.hasSigmas ? UI.sigma(table.orbit.sigmaPeriod) : "—"
                numeric: true
                secondary: true
                font.pixelSize: 11
                Layout.preferredWidth: 96
                horizontalAlignment: Text.AlignRight
                HoverHandler { id: periodSigmaHover }
                ToolTip.visible: periodSigmaHover.hovered
                ToolTip.delay: 500
                ToolTip.text: qsTr("Period uncertainty in days")
            }
        }
    }
    Rectangle {
        Layout.fillWidth: true
        implicitHeight: 1
        color: AppTheme.lineSoft
        visible: table.orbit.eccentricity < 1 && Number.isFinite(table.orbit.orbitalPeriod) && table.orbit.orbitalPeriod > 0
    }
    CustomTextLink { text: qsTr("All parameters and vectors"); alwaysUnderline: true; implicitHeight: 16; Layout.topMargin: 12; Layout.bottomMargin: 16; onClicked: details.open() }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.lineSoft }
    SectionHeading { title: qsTr("Solution and dynamics"); Layout.fillWidth: true; Layout.topMargin: 16; Layout.bottomMargin: 4 }
    RowLayout {
        Layout.fillWidth: true
        Layout.preferredHeight: 34
        CustomText { text: qsTr("Parameter"); secondary: true; font.pixelSize: 10; Layout.fillWidth: true }
        CustomText { text: qsTr("Value"); secondary: true; font.pixelSize: 10; horizontalAlignment: Text.AlignRight }
    }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.lineSoft }
    RowLayout {
        objectName: "orbit-quality"
        Layout.fillWidth: true
        Layout.preferredHeight: 42
        CustomText {
            text: qsTr("Orbit quality U")
            font.pixelSize: 11
            Layout.fillWidth: true
            HoverHandler { id: qualityHover }
            ToolTip.visible: qualityHover.hovered
            ToolTip.delay: 500
            ToolTip.text: qsTr("Estimate from the orbit covariance; may differ from the MPC uncertainty parameter.")
        }
        CustomText {
            objectName: "orbit-quality-value"
            text: UI.number(table.solution.orbitQualityU, 1)
            numeric: true
            font.pixelSize: 11
        }
    }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.lineSoft }
    NonGravParameters {
        objectName: "nominal-non-grav"
        solution: table.solution
        Layout.fillWidth: true
        Layout.bottomMargin: 16
    }
    Rectangle { Layout.fillWidth: true; Layout.bottomMargin: 16; implicitHeight: 1; color: AppTheme.lineSoft }
    ColumnLayout {
        id: photometryTable
        objectName: "nominal-photometry"
        readonly property bool compact: width < 380
        readonly property real labelWidth: compact ? 78 : 92
        readonly property real sigmaWidth: compact ? 64 : 84
        readonly property real columnSpacing: compact ? 4 : 10
        Layout.fillWidth: true
        spacing: 0
        SectionHeading { title: qsTr("Photometry"); Layout.fillWidth: true; Layout.bottomMargin: 8 }
        RowLayout {
            Layout.fillWidth: true
            Layout.preferredHeight: 34
            spacing: photometryTable.columnSpacing
            CustomText { text: qsTr("Parameter"); secondary: true; font.pixelSize: 10; Layout.preferredWidth: photometryTable.labelWidth }
            CustomText { text: qsTr("Value"); secondary: true; font.pixelSize: 10; Layout.fillWidth: true; horizontalAlignment: Text.AlignRight }
            CustomText { text: qsTr("σ phot."); secondary: true; font.pixelSize: 10; Layout.preferredWidth: photometryTable.sigmaWidth; horizontalAlignment: Text.AlignRight }
            CustomText { text: qsTr("σ orbit"); secondary: true; font.pixelSize: 10; Layout.preferredWidth: photometryTable.sigmaWidth; horizontalAlignment: Text.AlignRight }
        }
        Repeater {
            model: table.solution.cometPhotometry ? [
                {label: "M1", value: table.solution.cometTotalMagnitude, sigma: table.solution.cometTotalMagnitudePhotometricSigma, orbitSigma: table.solution.cometTotalMagnitudeSigma, digits: 9, unit: "mag"},
                {label: "M2", value: table.solution.cometNuclearMagnitude, sigma: table.solution.cometNuclearMagnitudePhotometricSigma, orbitSigma: table.solution.cometNuclearMagnitudeSigma, digits: 9, unit: "mag"},
                {label: "K1", value: table.solution.cometK1, digits: 2},
                {label: "K2", value: table.solution.cometK2, digits: 2},
                {label: "PHCOF", value: table.solution.cometNuclearPhaseCoefficient, digits: 3, unit: "mag/deg"}
            ] : [
                {label: "H", value: table.solution.absoluteMagnitude, sigma: table.solution.absoluteMagnitudePhotometricSigma, orbitSigma: table.solution.absoluteMagnitudeSigma, digits: 3, unit: "mag"},
                {label: "G", value: table.solution.asteroidSlopeG, digits: 3}
            ]
            Item {
                id: photometryRow
                required property var modelData
                objectName: "nominal-photometry-" + modelData.label
                Layout.fillWidth: true
                implicitHeight: 42
                Rectangle { width: parent.width; height: 1; color: AppTheme.lineSoft }
                RowLayout {
                    anchors.fill: parent
                    spacing: photometryTable.columnSpacing
                    RowLayout {
                        Layout.preferredWidth: photometryTable.labelWidth
                        spacing: 4
                        CustomText {
                            text: photometryRow.modelData.label
                            secondary: true
                            font.pixelSize: 11
                            HoverHandler { id: photometryHover }
                            ToolTip.visible: (photometryRow.modelData.label === "H" || photometryRow.modelData.label === "M1" || photometryRow.modelData.label === "M2") && photometryHover.hovered
                            ToolTip.delay: 500
                            ToolTip.text: qsTr("Photometric σ assumes independent magnitude errors and a fixed model and geometry. Orbital σ propagates orbital covariance.")
                        }
                        CustomText { text: photometryRow.modelData.unit || ""; secondary: true; font.pixelSize: 10 }
                        Item { Layout.fillWidth: true }
                    }
                    CustomText {
                        objectName: "nominal-photometry-value"
                        text: UI.number(photometryRow.modelData.value, photometryRow.modelData.digits)
                        numeric: true
                        font.pixelSize: 11
                        Layout.fillWidth: true
                        horizontalAlignment: Text.AlignRight
                        elide: Text.ElideLeft
                    }
                    CustomText {
                        text: UI.sigma(photometryRow.modelData.sigma)
                        numeric: true
                        secondary: true
                        font.pixelSize: photometryTable.compact ? 10 : 11
                        Layout.preferredWidth: photometryTable.sigmaWidth
                        horizontalAlignment: Text.AlignRight
                        elide: Text.ElideLeft
                    }
                    CustomText {
                        text: UI.sigma(photometryRow.modelData.orbitSigma)
                        numeric: true
                        secondary: true
                        font.pixelSize: photometryTable.compact ? 10 : 11
                        Layout.preferredWidth: photometryTable.sigmaWidth
                        horizontalAlignment: Text.AlignRight
                        elide: Text.ElideLeft
                    }
                }
            }
        }
        Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.lineSoft }
    }
    Rectangle {
        Layout.fillWidth: true
        Layout.topMargin: 16
        implicitHeight: entryContent.implicitHeight + 24
        color: UI.mix(AppTheme.amber, AppTheme.surface, 0.08)
        visible: table.solution.nominalEarthEntry.utc !== undefined
        Rectangle { width: 3; height: parent.height; color: AppTheme.amber }
        ColumnLayout {
            id: entryContent
            anchors.fill: parent
            anchors.leftMargin: 15
            anchors.rightMargin: 12
            anchors.topMargin: 12
            anchors.bottomMargin: 12
            spacing: 6
            RowLayout {
                Layout.fillWidth: true
                CustomText { text: qsTr("Nominal Earth entry"); color: AppTheme.amber; font.pixelSize: 12; Layout.fillWidth: true }
                CustomText { text: qsTr("%1 km").arg(table.solution.nominalEarthEntry.altitudeKm); secondary: true; font.pixelSize: 10 }
            }
            Repeater {
                model: [
                    { label: qsTr("Time, UTC"), value: table.solution.nominalEarthEntry.utc || "" },
                    { label: qsTr("Latitude"), value: table.solution.nominalEarthEntry.latitude || "" },
                    { label: qsTr("Longitude"), value: table.solution.nominalEarthEntry.longitude || "" }
                ]
                RowLayout {
                    required property var modelData
                    Layout.fillWidth: true
                    CustomText { text: modelData.label; secondary: true; font.pixelSize: 11; Layout.fillWidth: true }
                    CustomText { text: modelData.value; numeric: true; font.pixelSize: 11 }
                }
            }
        }
    }
    CustomText {
        text: qsTr("Nominal entry check failed: %1").arg(table.solution.nominalEarthEntryError)
        visible: table.solution.nominalEarthEntryError.length > 0
        Layout.fillWidth: true
        Layout.topMargin: 12
        wrapMode: Text.Wrap
        color: AppTheme.amber
        font.pixelSize: 11
    }
    MoidSummary { solution: table.solution; Layout.fillWidth: true; Layout.topMargin: 17 }
    OrbitDetailsDialog { id: details; solution: table.solution; elementRows: table.elementRows; equinoctial: table.equinoctial }
}
