import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import TrajectaUI

Rectangle {
    id: panel
    required property RelativeMotionModel model
    signal captureSelected(int index)
    function formatTimeUncertainty(seconds) {
        if (seconds === undefined) return qsTr("\u00b1 n/a")
        if (seconds < 1) return "\u00b1 " + seconds.toFixed(3) + " s"
        if (seconds < 60) return "\u00b1 " + seconds.toFixed(1) + " s"
        if (seconds < 3600) return "\u00b1 " + (seconds / 60).toFixed(1) + " min"
        if (seconds < 86400) return "\u00b1 " + (seconds / 3600).toFixed(1) + " h"
        return "\u00b1 " + (seconds / 86400).toFixed(2) + " d"
    }
    implicitWidth: 350
    color: AppTheme.surface
    Rectangle { anchors.left: parent.left; width: 1; height: parent.height; color: AppTheme.line }
    Flickable {
        id: panelScroll
        objectName: "relative-panel-scroll"
        anchors.fill: parent
        anchors.margins: 18
        anchors.rightMargin: 36
        contentHeight: content.implicitHeight
        boundsBehavior: Flickable.StopAtBounds
        clip: true
        ScrollBar.vertical: CustomScrollBar {
            id: panelScrollBar
            objectName: "relative-panel-scrollbar"
            parent: panel
            x: panel.width - width - 10
            y: panelScroll.y
            height: panelScroll.height
            width: 8
            policy: panelScroll.contentHeight > panelScroll.height + 1 ? ScrollBar.AlwaysOn : ScrollBar.AlwaysOff
            active: true
            contentItem: Rectangle { implicitWidth: 6; radius: 3; color: panelScrollBar.pressed ? AppTheme.ink : AppTheme.muted }
        }
        ColumnLayout {
            id: content
            width: parent.width
            spacing: 12
            CustomText { text: qsTr("Relative motion"); font.pixelSize: 15; font.weight: Font.DemiBold }
            CustomText { text: qsTr("Center"); secondary: true; font.pixelSize: 11 }
            CustomComboBox {
                id: planet
                objectName: "relative-planet"
                Layout.fillWidth: true
                model: ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
                currentIndex: panel.model.relativePlanetIndex
                onActivated: panel.model.relativePlanetIndex = currentIndex
                enabled: !panel.model.relativeBusy
                Accessible.name: qsTr("Central planet")
            }
            CustomText { text: qsTr("Time scale"); secondary: true; font.pixelSize: 11 }
            CustomComboBox {
                objectName: "relative-time-scale"
                Layout.fillWidth: true
                model: ["UTC", "TDB"]
                currentIndex: panel.model.relativeTimeScale === "TDB" ? 1 : 0
                onActivated: panel.model.relativeTimeScale = currentIndex === 1 ? "TDB" : "UTC"
                enabled: !panel.model.relativeBusy
                Accessible.name: qsTr("Time scale")
            }
            CustomText { text: qsTr("Start, %1").arg(panel.model.relativeTimeScale); secondary: true; font.pixelSize: 11 }
            CustomTextField { id: start; objectName: "relative-start"; Layout.fillWidth: true; text: panel.model.relativeStartDate; onTextEdited: panel.model.relativeStartDate = text; placeholderText: "YYYY-MM-DD"; enabled: !panel.model.relativeBusy }
            CustomText { text: qsTr("End, %1").arg(panel.model.relativeTimeScale); secondary: true; font.pixelSize: 11 }
            CustomTextField { id: end; objectName: "relative-end"; Layout.fillWidth: true; text: panel.model.relativeEndDate; onTextEdited: panel.model.relativeEndDate = text; placeholderText: "YYYY-MM-DD"; enabled: !panel.model.relativeBusy }
            CustomText { text: qsTr("Sampling step, days"); secondary: true; font.pixelSize: 11 }
            CustomNumberField { id: step; objectName: "relative-step"; Layout.fillWidth: true; value: panel.model.relativeStepDays; onValueEdited: number => panel.model.relativeStepDays = number; minimumValue: 0; exclusiveMinimum: true; enabled: !panel.model.relativeBusy }
            CustomText { text: qsTr("Reference frame"); secondary: true; font.pixelSize: 11 }
            CustomComboBox {
                id: frame
                objectName: "relative-frame"
                Layout.fillWidth: true
                model: [qsTr("Sun–planet rotating"), qsTr("Planet-centered inertial")]
                currentIndex: panel.model.rotatingFrame ? 0 : 1
                onActivated: panel.model.rotatingFrame = currentIndex === 0
            }
            CustomButton {
                objectName: "relative-calculate"
                Layout.fillWidth: true
                text: panel.model.relativeBusy ? qsTr("Cancel") : qsTr("Calculate trajectory")
                iconName: panel.model.relativeBusy ? "x" : "play"
                primary: !panel.model.relativeBusy
                onClicked: {
                    if (panel.model.relativeBusy) panel.model.cancelRelative()
                    else panel.model.calculateRelative(start.text, end.text, planet.currentIndex, step.value)
                }
            }
            CustomText { Layout.fillWidth: true; visible: panel.model.relativeError.length > 0; text: panel.model.relativeError; color: AppTheme.red; wrapMode: Text.Wrap; font.pixelSize: 11 }
            Rectangle { Layout.fillWidth: true; height: 1; color: AppTheme.line; Layout.topMargin: 8 }
            CustomText { text: qsTr("At selected epoch"); font.weight: Font.DemiBold }
            KeyValueList {
                Layout.fillWidth: true
                spacing: 10
                rows: panel.model.relativeAvailable ? [
					[qsTr("Center"), panel.model.relativeCenter],
                    [qsTr("Distance"), panel.model.relativeDistanceKm.toFixed(0) + " km"],
                    [qsTr("Hill radius"), (panel.model.relativeHillRadius * 149597870.7).toFixed(0) + " km"],
                    [qsTr("Two-body energy"), panel.model.relativeEnergy.toFixed(4) + " km²/s²"]
                ] : []
            }
            CustomText { Layout.fillWidth: true; visible: panel.model.relativeAvailable; text: qsTr("Values use the nearest trajectory sample."); secondary: true; wrapMode: Text.Wrap; font.pixelSize: 10 }
            Rectangle { Layout.fillWidth: true; visible: panel.model.relativeAvailable; height: 1; color: AppTheme.line; Layout.topMargin: 8 }
            CustomText { visible: panel.model.relativeAvailable; text: qsTr("Temporary capture"); font.weight: Font.DemiBold }
            CustomText {
                visible: panel.model.relativeAvailable && panel.model.relativeCaptures.length === 0
                Layout.fillWidth: true
                text: qsTr("No intervals found")
                secondary: true
                font.pixelSize: 11
            }
            Repeater {
                model: panel.model.relativeCaptures
                ItemDelegate {
                    id: captureEvent
                    required property var modelData
                    required property int index
                    objectName: "relative-capture-event-" + index
                    Layout.fillWidth: true
                    padding: 6
                    hoverEnabled: true
                    InteractionCursor {}
                    Accessible.name: qsTr("Go to start of event %1").arg(index + 1)
                    onClicked: panel.captureSelected(index)
                    background: Rectangle {
                        radius: 4
                        color: captureEvent.hovered ? AppTheme.nav : "transparent"
                        border.color: captureEvent.activeFocus ? AppTheme.accent : "transparent"
                    }
                    contentItem: ColumnLayout {
                    spacing: 8
                    CustomText {
                        text: qsTr("Event %1").arg(captureEvent.index + 1)
                        color: AppTheme.accent
                        font.weight: Font.DemiBold
                        font.pixelSize: 11
                    }
                    ColumnLayout {
                        id: captureRows
                        Layout.fillWidth: true
                        spacing: 8
                        property var rows: [
                            [qsTr("Start, %1").arg(panel.model.relativeResultTimeScale), (captureEvent.modelData.startClipped ? "\u2264 " : "") + captureEvent.modelData.start.slice(0, 19) + (captureEvent.modelData.startClipped ? "" : " " + panel.formatTimeUncertainty(captureEvent.modelData.start3SigmaSeconds))],
                            [qsTr("End, %1").arg(panel.model.relativeResultTimeScale), (captureEvent.modelData.endClipped ? "\u2265 " : "") + captureEvent.modelData.end.slice(0, 19) + (captureEvent.modelData.endClipped ? "" : " " + panel.formatTimeUncertainty(captureEvent.modelData.end3SigmaSeconds))],
                            [qsTr("Duration"), (captureEvent.modelData.startClipped || captureEvent.modelData.endClipped ? "\u2265 " : "") + captureEvent.modelData.durationDays.toFixed(2) + " d"],
                            [qsTr("Inertial net"), captureEvent.modelData.inertialRevolutionsValid ? captureEvent.modelData.inertialRevolutions.toFixed(2) + " " + qsTr("turns") : qsTr("Unavailable")],
                            [qsTr("Inertial travel"), captureEvent.modelData.inertialRevolutionsValid ? captureEvent.modelData.inertialAngularTravel.toFixed(2) + " " + qsTr("turns") : qsTr("Unavailable")],
                            [qsTr("Rotating net"), captureEvent.modelData.revolutionsValid ? captureEvent.modelData.revolutions.toFixed(2) + " " + qsTr("turns") : qsTr("Unavailable")],
                            [qsTr("Rotating travel"), captureEvent.modelData.revolutionsValid ? captureEvent.modelData.angularTravel.toFixed(2) + " " + qsTr("turns") : qsTr("Unavailable")],
                            [qsTr("Coverage"), captureEvent.modelData.startClipped || captureEvent.modelData.endClipped ? qsTr("Partial interval") : qsTr("Complete interval")]
                        ]
                        Repeater {
                            model: captureRows.rows
                            RowLayout {
                                id: captureRow
                                required property var modelData
                                required property int index
                                Layout.fillWidth: true
                                spacing: 6
                                CustomText { text: modelData[0]; secondary: true; font.pixelSize: 11; Layout.preferredWidth: 90 }
                                CustomText { objectName: "relative-capture-value-" + captureRow.index; text: captureRow.modelData[1]; numeric: true; font.pixelSize: 10; font.weight: Font.DemiBold; fontSizeMode: Text.HorizontalFit; minimumPixelSize: 8; horizontalAlignment: Text.AlignLeft; Layout.fillWidth: true; Layout.minimumWidth: 0 }
                            }
                        }
                    }
                    }
                }
            }
        }
    }
}
