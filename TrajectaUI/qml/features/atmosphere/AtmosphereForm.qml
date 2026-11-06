pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import QtQuick.Dialogs
import QtQuick.Controls.Basic
import TrajectaUI

ColumnLayout {
    id: form
    objectName: "atmosphere-form"
    required property AtmosphereViewModel atmosphere
    spacing: 20
    component Field: ColumnLayout {
        id: field
        property string label
        Layout.fillWidth: true
        Layout.preferredWidth: 1
        spacing: 6
        CustomText { text: field.label; secondary: true; font.pixelSize: 11; Layout.fillWidth: true; wrapMode: Text.Wrap }
    }
    component Heading: RowLayout {
        id: heading
        property string label
        Layout.fillWidth: true
        Layout.topMargin: 10
        spacing: 12
        CustomText { text: heading.label; font.pixelSize: 11; font.weight: Font.DemiBold }
        Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
    }
    component Fields: GridLayout {
        Layout.fillWidth: true
        columns: width < 460 ? 1 : 2
        uniformCellWidths: true
        columnSpacing: 20
        rowSpacing: 18
    }
    RowLayout {
        Layout.fillWidth: true
        spacing: 12
        CustomText { text: qsTr("Atmospheric MC"); font.pixelSize: 18; font.weight: Font.DemiBold }
        CustomText { text: "StrewnField"; color: AppTheme.accent; font.pixelSize: 10 }
        Item { Layout.fillWidth: true }
    }
    RowLayout {
        Layout.fillWidth: true
        Layout.topMargin: -10
        spacing: 9
        CustomIcon { name: "file-text"; color: AppTheme.muted; size: 16 }
        CustomText { text: qsTr("Earth entries at 100 km / fixed entry"); secondary: true; font.pixelSize: 11 }
    }
    Heading { label: qsTr("Initial state") }
    Fields {
        Field {
            label: qsTr("Source")
            CustomComboBox {
                objectName: "atmosphere-source-mode"
                Layout.fillWidth: true
                model: [qsTr("Entry archive"), qsTr("Fixed entry")]
                currentIndex: form.atmosphere.sourceMode
                Accessible.name: qsTr("Atmospheric entry source")
                onActivated: form.atmosphere.sourceMode = currentIndex
            }
        }
        Field {
            visible: form.atmosphere.sourceMode === 0
            label: qsTr("Orbital run")
            CustomComboBox {
                objectName: "atmosphere-orbital-run"
                Layout.fillWidth: true
                model: [{label: qsTr("Select a run")}].concat(form.atmosphere.sources)
                textRole: "label"
                currentIndex: {
                    for (let i = 0; i < form.atmosphere.sources.length; ++i)
                        if (form.atmosphere.sources[i].path === form.atmosphere.archivePath) return i + 1
                    return 0
                }
                enabled: form.atmosphere.sources.length > 0
                Accessible.name: qsTr("Orbital MC run at 100 km")
                onActivated: if (currentIndex > 0) form.atmosphere.selectArchive(form.atmosphere.sources[currentIndex - 1].url)
            }
        }
        AtmosphereNumberField {
            visible: form.atmosphere.sourceMode === 1
            atmosphere: form.atmosphere
            settingKey: "epoch"
            label: qsTr("Entry epoch, UTC")
            placeholder: "YYYY-MM-DD HH:MM:SS"
        }
        Field {
            visible: form.atmosphere.sourceMode === 0
            Layout.columnSpan: parent.columns
            label: qsTr("Entry archive")
            RowLayout {
                Layout.fillWidth: true
                CustomTextField {
                    objectName: "atmosphere-archive"
                    Layout.fillWidth: true
                    readOnly: true
                    text: form.atmosphere.archivePath
                    placeholderText: qsTr("No archive selected")
                    Accessible.name: qsTr("Selected entry archive")
                }
                CustomButton {
                    objectName: "atmosphere-open-archive"
                    iconName: "folder-open"
                    accessibleLabel: qsTr("Open entry archive")
                    onClicked: archiveDialog.open()
                    ToolTip.visible: hovered
                    ToolTip.text: qsTr("Open entry archive")
                }
            }
        }
        Repeater {
            model: [
                {key: "latitude", label: qsTr("Latitude WGS84, degrees")},
                {key: "longitude", label: qsTr("Longitude WGS84, degrees")},
                {key: "speed", label: qsTr("Speed ECEF, km/s")},
                {key: "heading", label: qsTr("Heading from north, degrees")},
                {key: "slope", label: qsTr("Slope from horizon, degrees")}
            ]
            AtmosphereNumberField {
                required property var modelData
                visible: form.atmosphere.sourceMode === 1
                atmosphere: form.atmosphere
                settingKey: modelData.key
                label: modelData.label
            }
        }
        Field {
            label: qsTr("Boundary")
            CustomText { text: qsTr("Earth / WGS84 / 100 km"); font.pixelSize: 12; Layout.fillWidth: true; wrapMode: Text.Wrap }
        }
        Field {
            visible: form.atmosphere.sourceMode === 0 && form.atmosphere.events.length > 0
            Layout.columnSpan: parent.columns
            label: qsTr("Entry event, UTC")
            CustomComboBox {
                objectName: "atmosphere-event"
                Layout.fillWidth: true
                model: form.atmosphere.events
                currentIndex: form.atmosphere.selectedEvent
                displayText: currentIndex < 0 ? qsTr("Select one event") : currentText
                onActivated: form.atmosphere.selectedEvent = currentIndex
                Accessible.name: qsTr("Atmospheric entry event")
            }
        }
    }
    CustomText {
        objectName: "atmosphere-source-error"
        text: form.atmosphere.sourceError
        visible: text.length > 0
        color: AppTheme.amber
        font.pixelSize: 11
        wrapMode: Text.Wrap
        elide: Text.ElideNone
        Layout.fillWidth: true
    }
    CustomText {
        text: form.atmosphere.sourceSummary
        visible: text.length > 0
        secondary: true
        font.pixelSize: 11
        Layout.fillWidth: true
        wrapMode: Text.Wrap
        elide: Text.ElideNone
    }
    Heading { label: qsTr("Body and fragmentation") }
    Fields {
        Field {
            label: qsTr("Material preset")
            CustomComboBox {
                objectName: "atmosphere-preset"
                Layout.fillWidth: true
                model: [qsTr("Ordinary chondrite"), qsTr("Carbonaceous chondrite"), qsTr("Iron"), qsTr("Custom")]
                currentIndex: form.atmosphere.preset
                Accessible.name: qsTr("Material preset")
                onActivated: {
                    if (currentIndex < 3) form.atmosphere.selectPreset(currentIndex)
                    currentIndex = Qt.binding(() => form.atmosphere.preset)
                }
            }
        }
        Field {
            label: qsTr("Initial size")
            CustomComboBox {
                objectName: "atmosphere-size-mode"
                Layout.fillWidth: true
                model: [qsTr("H and albedo"), qsTr("Mass range"), qsTr("Diameter range")]
                currentIndex: form.atmosphere.sizeMode
                Accessible.name: qsTr("Initial size distribution")
                onActivated: form.atmosphere.sizeMode = currentIndex
            }
        }
        AtmosphereNumberField { atmosphere: form.atmosphere; settingKey: "magnitude"; label: qsTr("Absolute magnitude H"); visible: form.atmosphere.sizeMode === 0 }
        Field {
            visible: form.atmosphere.sizeMode === 0
            label: qsTr("Albedo, min / mode / max")
            RowLayout {
                Layout.fillWidth: true
                spacing: 8
                AtmosphereNumberField { atmosphere: form.atmosphere; settingKey: "albedoMin"; placeholder: qsTr("Minimum albedo") }
                AtmosphereNumberField { atmosphere: form.atmosphere; settingKey: "albedoMode"; placeholder: qsTr("Modal albedo") }
                AtmosphereNumberField { atmosphere: form.atmosphere; settingKey: "albedoMax"; placeholder: qsTr("Maximum albedo") }
            }
        }
        AtmosphereRangeField { atmosphere: form.atmosphere; settingKey: "mass"; label: qsTr("Initial mass, kg"); visible: form.atmosphere.sizeMode === 1 }
        AtmosphereRangeField { atmosphere: form.atmosphere; settingKey: "diameter"; label: qsTr("Initial diameter, m"); visible: form.atmosphere.sizeMode === 2 }
        Repeater {
            model: [
                {key: "density", label: qsTr("Density, kg/m3")},
                {key: "strength", label: qsTr("Strength, MPa (q = rho v^2)")},
                {key: "ablation", label: qsTr("Ablation, 1e-8 kg/J")},
                {key: "dust", label: qsTr("Dust fraction")},
                {key: "largest", label: qsTr("Largest fragment fraction")},
                {key: "beta", label: qsTr("Fragment mass exponent beta")}
            ]
            AtmosphereRangeField {
                required property var modelData
                atmosphere: form.atmosphere
                settingKey: modelData.key
                label: modelData.label
            }
        }
        AtmosphereNumberField { atmosphere: form.atmosphere; settingKey: "alpha"; label: qsTr("Weibull alpha") }
        AtmosphereNumberField { atmosphere: form.atmosphere; settingKey: "minMass"; label: qsTr("Minimum fragment, g") }
        AtmosphereNumberField { atmosphere: form.atmosphere; settingKey: "representatives"; label: qsTr("Representatives per fragment bin") }
    }
    Heading { label: qsTr("Environment") }
    Fields {
        Field {
            label: qsTr("Atmosphere and terrain")
            CustomComboBox {
                objectName: "atmosphere-environment-mode"
                Layout.fillWidth: true
                model: [qsTr("Automatic GFS + DEM"), qsTr("Saved environment"), qsTr("Idealized atmosphere + flat ground")]
                currentIndex: form.atmosphere.environmentMode
                Accessible.name: qsTr("Environment model")
                onActivated: form.atmosphere.environmentMode = currentIndex
            }
        }
        CustomCheckBox {
            objectName: "atmosphere-offline"
            Layout.fillWidth: true
            Layout.preferredWidth: 1
            Layout.alignment: Qt.AlignBottom
            text: qsTr("Offline / cached data only")
            visible: form.atmosphere.environmentMode === 0
            checked: form.atmosphere.offline
            onToggled: form.atmosphere.offline = checked
        }
        Field {
            visible: form.atmosphere.environmentMode === 1
            label: qsTr("Environment profile")
            RowLayout {
                Layout.fillWidth: true
                CustomTextField { Layout.fillWidth: true; readOnly: true; text: form.atmosphere.environmentPath; placeholderText: qsTr("No profile selected") }
                CustomButton {
                    iconName: "folder-open"
                    accessibleLabel: qsTr("Open environment profile")
                    onClicked: environmentDialog.open()
                    ToolTip.visible: hovered
                    ToolTip.text: qsTr("Open environment profile")
                }
            }
        }
        AtmosphereNumberField { atmosphere: form.atmosphere; settingKey: "groundHeight"; label: qsTr("Flat ground height WGS84, m"); visible: form.atmosphere.environmentMode === 2 }
        AtmosphereRangeField { atmosphere: form.atmosphere; settingKey: "wind"; label: qsTr("Wind multiplier") }
        AtmosphereRangeField { atmosphere: form.atmosphere; settingKey: "drag"; label: qsTr("Drag multiplier") }
    }
    Heading { label: qsTr("Sampling and convergence") }
    Fields {
        Repeater {
            model: [
                {key: "runs", label: qsTr("Maximum realizations")},
                {key: "seed", label: qsTr("Seed")},
                {key: "workers", label: qsTr("Workers")},
                {key: "timeout", label: qsTr("Realization timeout, s")}
            ]
            AtmosphereNumberField { required property var modelData; atmosphere: form.atmosphere; settingKey: modelData.key; label: modelData.label }
        }
    }
    CustomCheckBox { objectName: "atmosphere-convergence"; text: qsTr("Stop on JS convergence"); checked: form.atmosphere.convergence; onToggled: form.atmosphere.convergence = checked }
    Fields {
        visible: form.atmosphere.convergence
        AtmosphereNumberField { atmosphere: form.atmosphere; settingKey: "threshold"; label: qsTr("JS divergence threshold") }
        AtmosphereNumberField { atmosphere: form.atmosphere; settingKey: "batch"; label: qsTr("Check interval, realizations") }
        AtmosphereNumberField { atmosphere: form.atmosphere; settingKey: "checks"; label: qsTr("Consecutive stable checks") }
    }
    Heading { label: qsTr("Output") }
    Fields {
        Field {
            label: qsTr("Map smoothing")
            CustomComboBox {
                objectName: "atmosphere-map-method"
                Layout.fillWidth: true
                model: [qsTr("Adaptive KDE"), qsTr("Fixed KDE")]
                currentIndex: form.atmosphere.mapMethod
                Accessible.name: qsTr("Map smoothing method")
                onActivated: form.atmosphere.mapMethod = currentIndex
            }
        }
        AtmosphereNumberField { atmosphere: form.atmosphere; settingKey: "cell"; label: qsTr("Map cell size, m") }
        AtmosphereNumberField { atmosphere: form.atmosphere; settingKey: "sigma"; label: qsTr("Kernel sigma, m"); visible: form.atmosphere.mapMethod === 1 }
    }
    RowLayout {
        spacing: 12
        CustomIcon { name: "file-text"; color: AppTheme.blue; size: 18 }
        ColumnLayout {
            Layout.fillWidth: true
            spacing: 7
            CustomText { text: qsTr("Numerical archive + HTML report"); font.pixelSize: 11; font.weight: Font.DemiBold; Layout.fillWidth: true; wrapMode: Text.Wrap }
            CustomText { text: qsTr("Not generated"); secondary: true; font.pixelSize: 10 }
        }
    }
    CustomText { objectName: "atmosphere-form-error"; text: form.atmosphere.formError; visible: text.length > 0; color: AppTheme.red; font.pixelSize: 11; Layout.fillWidth: true; wrapMode: Text.Wrap; elide: Text.ElideNone }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
    CustomText { text: form.atmosphere.error; visible: text.length > 0; color: AppTheme.red; font.pixelSize: 11; Layout.fillWidth: true; wrapMode: Text.Wrap; elide: Text.ElideNone }
    FileDialog {
        id: archiveDialog
        title: qsTr("Select orbital MC result.json")
        nameFilters: [qsTr("MC archive (result.json)"), qsTr("JSON files (*.json)")]
        onAccepted: form.atmosphere.selectArchive(selectedFile)
    }
    FileDialog {
        id: environmentDialog
        title: qsTr("Select environment profile")
        nameFilters: [qsTr("Environment files (*.json)")]
        onAccepted: form.atmosphere.selectEnvironment(selectedFile)
    }
}
