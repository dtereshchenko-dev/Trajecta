pragma ComponentBehavior: Bound
import QtCore
import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Dialogs
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

BasePage {
    id: page
    objectName: "ephemerides-page"
    required property EphemerisViewModel ephemerides
    readonly property real pageMargin: UI.pageMargin(width)
    signal importRequested()
    title: qsTr("Ephemerides")
    padding: 0
    header: Rectangle {
        implicitHeight: 48
        color: AppTheme.base
        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: page.pageMargin
            anchors.rightMargin: page.pageMargin
            spacing: 14
            CustomText { text: page.title; font.pixelSize: 14; font.weight: Font.DemiBold; Layout.fillWidth: true }
            CustomText {
                objectName: "ephemeris-planned-epochs"
                visible: page.ephemerides.available
                text: page.ephemerides.plannedEpochs === 1 ? qsTr("1 epoch") : qsTr("%1 epochs").arg(page.ephemerides.plannedEpochs)
                secondary: true
                font.pixelSize: 11
            }
            CustomButton {
                objectName: "ephemeris-settings-button"
                iconName: "sliders-horizontal"
                flat: true
                accessibleLabel: qsTr("Time scale and light time")
                onClicked: settings.open()
            }
        }
        Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
    }
    contentItem: Item {
        ColumnLayout {
            anchors.fill: parent
            anchors.margins: page.pageMargin
            spacing: 16
            RowLayout {
                Layout.fillWidth: true
                spacing: 12
                CustomSegmentedControl {
                    objectName: "ephemeris-layout"
                    options: [qsTr("Spherical"), qsTr("Cartesian")]
                    currentIndex: page.ephemerides.layout
                    onActivated: index => page.ephemerides.layout = index
                }
                Item { Layout.fillWidth: true }
                CustomSegmentedControl {
                    objectName: "ephemeris-center"
                    options: [qsTr("Helio"), qsTr("Geo"), qsTr("Topo")]
                    currentIndex: page.ephemerides.center
                    onActivated: index => page.ephemerides.center = index
                }
                CustomText { visible: page.ephemerides.center === 2; text: qsTr("Station"); secondary: true; font.pixelSize: 11 }
                CustomTextField {
                    objectName: "ephemeris-station"
                    visible: page.ephemerides.center === 2
                    implicitWidth: 64
                    font.family: UI.monoFont
                    font.pixelSize: 12
                    maximumLength: 3
                    text: page.ephemerides.station
                    color: text.trim().length ? AppTheme.ink : AppTheme.red
                    Accessible.name: qsTr("Station code")
                    onTextEdited: page.ephemerides.station = text
                }
            }
            RowLayout {
                Layout.fillWidth: true
                spacing: 14
                ColumnLayout {
                    spacing: 5
                    CustomText { text: qsTr("Start · %1").arg(page.ephemerides.scaleName); secondary: true; font.pixelSize: 10 }
                    CustomTextField {
                        objectName: "ephemeris-start"
                        implicitWidth: 124
                        font.family: UI.monoFont
                        font.pixelSize: 12
                        placeholderText: "yyyy-mm-dd"
                        text: page.ephemerides.startDate
                        color: page.ephemerides.startValid ? AppTheme.ink : AppTheme.red
                        Accessible.name: qsTr("Start date")
                        onTextEdited: page.ephemerides.startDate = text
                    }
                }
                ColumnLayout {
                    spacing: 5
                    CustomText { text: qsTr("End · %1").arg(page.ephemerides.scaleName); secondary: true; font.pixelSize: 10 }
                    CustomTextField {
                        objectName: "ephemeris-end"
                        implicitWidth: 124
                        font.family: UI.monoFont
                        font.pixelSize: 12
                        placeholderText: "yyyy-mm-dd"
                        text: page.ephemerides.endDate
                        color: page.ephemerides.endValid ? AppTheme.ink : AppTheme.red
                        Accessible.name: qsTr("End date")
                        onTextEdited: page.ephemerides.endDate = text
                    }
                }
                ColumnLayout {
                    spacing: 5
                    CustomText { text: qsTr("Step, days"); secondary: true; font.pixelSize: 10 }
                    CustomNumberField {
                        objectName: "ephemeris-step"
                        implicitWidth: 96
                        value: page.ephemerides.stepDays
                        minimumValue: 0
                        exclusiveMinimum: true
                        Accessible.name: qsTr("Step in days")
                        onValueEdited: number => page.ephemerides.stepDays = number
                    }
                }
                ColumnLayout {
                    spacing: 5
                    CustomText { text: qsTr("Plane"); secondary: true; font.pixelSize: 10 }
                    CustomComboBox {
                        objectName: "ephemeris-plane"
                        implicitWidth: 150
                        model: [qsTr("Equatorial"), qsTr("Ecliptic")]
                        currentIndex: page.ephemerides.plane
                        Accessible.name: qsTr("Reference plane")
                        onActivated: page.ephemerides.plane = currentIndex
                    }
                }
                Item { Layout.fillWidth: true }
                CustomButton {
                    objectName: "ephemeris-generate"
                    Layout.alignment: Qt.AlignBottom
                    text: page.ephemerides.busy ? qsTr("Cancel") : qsTr("Generate")
                    iconName: page.ephemerides.busy ? "square" : "play"
                    primary: !page.ephemerides.busy
                    enabled: page.ephemerides.busy || page.ephemerides.canGenerate
                    onClicked: page.ephemerides.busy ? page.ephemerides.cancel() : page.ephemerides.generate()
                }
            }
            CustomText {
                objectName: "ephemeris-form-error"
                visible: text.length > 0
                text: page.ephemerides.formError
                color: AppTheme.red
                font.pixelSize: 11
                Layout.fillWidth: true
            }
            RowLayout {
                Layout.fillWidth: true
                spacing: 9
                CustomText {
                    objectName: "ephemeris-summary"
                    text: page.ephemerides.busy ? qsTr("Computing ephemeris...") : page.ephemerides.summary
                    secondary: true
                    font.pixelSize: 11
                    Layout.fillWidth: true
                }
                CustomComboBox {
                    id: exportFormat
                    objectName: "ephemeris-export-format"
                    implicitWidth: 86
                    model: page.ephemerides.exportFormats
                    Accessible.name: qsTr("Export format")
                }
                CustomButton {
                    objectName: "ephemeris-export"
                    iconName: "download"
                    enabled: page.ephemerides.epochCount > 0 && !page.ephemerides.busy
                    accessibleLabel: qsTr("Save ephemerides")
                    onClicked: {
                        saveDialog.selectedFile = saveDialog.currentFolder.toString().replace(/\/$/, "") + "/"
                            + page.ephemerides.suggestedFileName(exportFormat.currentIndex)
                        saveDialog.open()
                    }
                }
            }
            Item {
                Layout.fillWidth: true
                Layout.fillHeight: true
                HorizontalHeaderView {
                    id: header
                    anchors.left: parent.left
                    anchors.right: parent.right
                    anchors.top: parent.top
                    anchors.rightMargin: verticalBar.visible ? verticalBar.width : 0
                    visible: table.rows > 0
                    syncView: table
                    clip: true
                    delegate: Rectangle {
                        id: headerCell
                        required property var model
                        required property int column
                        implicitHeight: 38
                        color: AppTheme.base
                        CustomText {
                            anchors.fill: parent
                            anchors.leftMargin: 8
                            anchors.rightMargin: 8
                            text: headerCell.model.display
                            color: AppTheme.muted
                            font.pixelSize: 10
                            horizontalAlignment: headerCell.column > 0 ? Text.AlignRight : Text.AlignLeft
                        }
                        Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
                    }
                }
                TableView {
                    id: table
                    objectName: "ephemeris-table"
                    anchors.left: parent.left
                    anchors.right: parent.right
                    anchors.top: header.bottom
                    anchors.bottom: parent.bottom
                    clip: true
                    model: page.ephemerides.table
                    boundsBehavior: Flickable.StopAtBounds
                    columnWidthProvider: function(column) {
                        const widths = page.ephemerides.layout === 1 ? [176, 128, 128, 128, 128, 128, 128]
                            : [176, 132, 124, 120, 100, 76, 108]
                        const total = widths.reduce((a, b) => a + b, 0)
                        return widths[column] + Math.max(0, width - total - 14) / widths.length
                    }
                    onWidthChanged: forceLayout()
                    rowHeightProvider: function(row) { return 34 }
                    delegate: Rectangle {
                        id: cell
                        required property var display
                        required property int column
                        implicitHeight: 34
                        color: AppTheme.base
                        CustomText {
                            anchors.fill: parent
                            anchors.leftMargin: 8
                            anchors.rightMargin: 8
                            text: cell.display === undefined || cell.display === null ? "—" : cell.display
                            numeric: true
                            font.pixelSize: 11
                            horizontalAlignment: cell.column > 0 ? Text.AlignRight : Text.AlignLeft
                        }
                        Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.lineSoft }
                    }
                    ScrollBar.horizontal: CustomScrollBar { id: horizontalBar; policy: ScrollBar.AsNeeded }
                    ScrollBar.vertical: CustomScrollBar { id: verticalBar; policy: ScrollBar.AsNeeded }
                    bottomMargin: horizontalBar.visible ? horizontalBar.height : 0
                    rightMargin: verticalBar.visible ? verticalBar.width : 0
                }
                ColumnLayout {
                    anchors.centerIn: parent
                    spacing: 16
                    visible: table.rows === 0
                    CustomIcon { Layout.alignment: Qt.AlignHCenter; name: "calendar-range"; size: 28; color: AppTheme.muted }
                    CustomText {
                        objectName: "ephemeris-empty"
                        Layout.alignment: Qt.AlignHCenter
                        text: !page.ephemerides.available ? qsTr("No solution")
                            : page.ephemerides.busy ? qsTr("Computing ephemeris...")
                            : qsTr("Choose the range and generate the ephemeris")
                    }
                    CustomButton {
                        objectName: "ephemeris-empty-action"
                        Layout.alignment: Qt.AlignHCenter
                        visible: !page.ephemerides.available
                        text: qsTr("Open observations")
                        iconName: "folder-open"
                        onClicked: page.importRequested()
                    }
                }
            }
            CustomText {
                objectName: "ephemeris-error"
                visible: text.length > 0
                text: page.ephemerides.error
                color: AppTheme.red
                wrapMode: Text.Wrap
                elide: Text.ElideNone
                font.pixelSize: 11
                Layout.fillWidth: true
            }
        }
    }
    Connections {
        target: page.ephemerides
        function onFormChanged() { table.forceLayout() }
    }
    FileDialog {
        id: saveDialog
        objectName: "ephemeris-save-dialog"
        title: qsTr("Save ephemerides")
        fileMode: FileDialog.SaveFile
        currentFolder: StandardPaths.writableLocation(StandardPaths.DocumentsLocation)
        nameFilters: [[qsTr("Text (*.txt)"), qsTr("JSON (*.json)")][exportFormat.currentIndex], qsTr("All files (*)")]
        onAccepted: page.ephemerides.save(selectedFile, exportFormat.currentIndex)
    }
    EphemerisSettingsDialog { id: settings; ephemerides: page.ephemerides }
}
