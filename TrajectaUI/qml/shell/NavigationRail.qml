pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI
import "../shared/UI.js" as UI

Rectangle {
    id: rail
    objectName: "navigation-rail"
    required property ApplicationSession session
    property int currentPage: 0
    property bool compact: false
    property bool collapsible: true
    signal toggleCollapsed()
    signal navigate(int page)
    signal importRequested()
    signal settingsRequested()
    signal dataSourcesRequested()
    color: AppTheme.nav
    implicitWidth: compact ? UI.navigationRailWidth : UI.navigationWidth
    ColumnLayout {
        anchors.fill: parent
        anchors.leftMargin: rail.compact ? 8 : 14
        anchors.rightMargin: rail.compact ? 8 : 14
        anchors.topMargin: 22
        anchors.bottomMargin: 12
        spacing: 0
        RowLayout {
            visible: !rail.compact
            Layout.leftMargin: 9
            Layout.bottomMargin: 28
            spacing: 11
            CustomIcon { name: "orbit"; color: AppTheme.accent; size: 31 }
            ColumnLayout {
                spacing: 0
                CustomText { text: "Trajecta"; font.pixelSize: 19; font.weight: 650 }
            }
        }
        CustomButton {
            objectName: "compact-open-observations"
            visible: rail.compact
            Layout.fillWidth: true
            Layout.preferredHeight: 43
            Layout.bottomMargin: 20
            flat: true
            accessibleLabel: qsTr("Open observations")
            enabled: !rail.session.busy
            contentItem: Item { CustomIcon { anchors.centerIn: parent; name: "orbit"; size: 31; color: AppTheme.accent } }
            onClicked: rail.importRequested()
        }
        Button {
            id: picker
            objectName: "open-observations-button"
            visible: !rail.compact
            Layout.fillWidth: true
            Layout.preferredHeight: 66
            Layout.bottomMargin: 26
            padding: 10
            hoverEnabled: true
            enabled: !rail.session.busy
            InteractionCursor {}
            contentItem: RowLayout {
                spacing: 10
                CustomIcon { name: "circle-dot"; color: AppTheme.accent }
                ColumnLayout {
                    Layout.fillWidth: true
                    spacing: 3
                    CustomText { text: rail.session.solution.available ? rail.session.solution.objectTitle : qsTr("Open observations"); font.pixelSize: 14; font.weight: Font.DemiBold; Layout.fillWidth: true }
                    CustomText { text: rail.session.solution.available ? rail.session.solution.designation : rail.session.inputName || qsTr("No file selected"); font.pixelSize: 10; numeric: true; secondary: true; Layout.fillWidth: true }
                }
            }
            background: Rectangle { radius: 5; color: picker.hovered ? AppTheme.hover : AppTheme.surfaceSecondary; border.color: picker.activeFocus ? AppTheme.accent : AppTheme.controlLine }
            Accessible.name: qsTr("Open observations")
            onClicked: rail.importRequested()
        }
        ColumnLayout {
            Layout.fillWidth: true
            spacing: rail.compact ? 9 : 5
            Repeater {
                model: [{label: qsTr("Solution"), icon: "scan-line"}, {label: qsTr("Observations"), icon: "table-2"}, {label: qsTr("Orbit 3D"), icon: "orbit"}, {label: qsTr("Ephemerides"), icon: "calendar-range"}, {label: qsTr("Risk Monitor"), icon: "activity"}, {label: qsTr("Close Approaches"), icon: "search"}]
                NavigationItem {
                    required property var modelData
                    required property int index
                    objectName: "navigation-page-" + index
                    Layout.fillWidth: true
                    text: modelData.label
                    iconName: modelData.icon
                    compact: rail.compact
                    checked: rail.currentPage === index
                    count: index === 1 && rail.session.solution.available ? String(rail.session.solution.total) : ""
                    onClicked: rail.navigate(index)
                }
            }
        }
        Item { Layout.fillHeight: true }
        NavigationItem { objectName: "settings-button"; compact: rail.compact; Layout.fillWidth: true; text: qsTr("Settings"); iconName: "sliders-horizontal"; onClicked: rail.settingsRequested() }
        NavigationItem { objectName: "data-sources-button"; compact: rail.compact; Layout.fillWidth: true; text: qsTr("Data sources"); iconName: "table-2"; warning: rail.session.dataSources.hasWarning; onClicked: rail.dataSourcesRequested() }
        NavigationItem {
            objectName: "toggle-navigation"
            visible: rail.collapsible
            compact: rail.compact
            Layout.fillWidth: true
            text: rail.compact ? qsTr("Expand sidebar") : qsTr("Collapse sidebar")
            iconName: rail.compact ? "panel-left-open" : "panel-left-close"
            onClicked: rail.toggleCollapsed()
        }
        Rectangle { visible: !rail.compact; Layout.fillWidth: true; Layout.topMargin: 16; Layout.bottomMargin: 16; implicitHeight: 1; color: AppTheme.line }
        RowLayout {
            visible: !rail.compact
            Layout.leftMargin: 4
            spacing: 6
            Rectangle { width: 6; height: 6; radius: 3; color: AppTheme.accent }
            CustomText { text: qsTr("Local session"); secondary: true; font.pixelSize: 10; Layout.fillWidth: true }
        }
    }
    Rectangle { anchors.right: parent.right; height: parent.height; width: 1; color: AppTheme.line }
}
