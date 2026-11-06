pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    id: panel
    required property AtmosphereViewModel atmosphere
    spacing: 20
    CustomText { text: qsTr("Source data"); font.pixelSize: 11; font.weight: Font.DemiBold }
    KeyValueList {
        Layout.fillWidth: true
        rows: [[qsTr("Source"), panel.atmosphere.sourceMode === 0 ? qsTr("Entry archive") : qsTr("Fixed entry")],
            [qsTr("Target"), qsTr("Earth")], [qsTr("Entry height"), "100 km"],
            [qsTr("Coordinates"), "WGS84"], [qsTr("Velocity frame"), "ECEF"], [qsTr("Date scale"), "UTC"],
            [qsTr("Size"), [qsTr("H and albedo"), qsTr("Mass range"), qsTr("Diameter range")][panel.atmosphere.sizeMode]],
            [qsTr("Environment"), ["GFS + DEM", qsTr("Saved profile"), qsTr("Idealized")][panel.atmosphere.environmentMode]],
            [qsTr("Probability"), qsTr("Conditional on entry")]]
    }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
    CustomText { text: qsTr("Route"); font.pixelSize: 11; font.weight: Font.DemiBold }
    Repeater {
        model: [qsTr("Entry states / fixed entry"), qsTr("Weather and terrain"), qsTr("Atmospheric realizations"), qsTr("Strewn field and HTML")]
        RowLayout {
            id: route
            required property string modelData
            required property int index
            Layout.fillWidth: true
            spacing: 16
            CustomText { text: "0" + (route.index + 1); numeric: true; color: AppTheme.accent; font.pixelSize: 10 }
            CustomText { text: route.modelData; secondary: true; font.pixelSize: 11; Layout.fillWidth: true; wrapMode: Text.Wrap; elide: Text.ElideNone }
        }
    }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
    CustomText { text: qsTr("Outcomes"); font.pixelSize: 11; font.weight: Font.DemiBold }
    CustomText { text: qsTr("Ground impact / Ablation / Incomplete / Failure"); secondary: true; font.pixelSize: 11; Layout.fillWidth: true; wrapMode: Text.Wrap; elide: Text.ElideNone }
}
