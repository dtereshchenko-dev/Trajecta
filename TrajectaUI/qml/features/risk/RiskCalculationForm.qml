pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: form
    required property RiskMonitorViewModel monitor
    spacing: 20
    component Field: ColumnLayout {
        id: field
        property string label
        spacing: 6
        Layout.fillWidth: true
        Layout.preferredWidth: 1
        CustomText { text: field.label; secondary: true; font.pixelSize: 11 }
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
    RowLayout {
        Layout.fillWidth: true
        spacing: 12
        CustomText { text: form.monitor.method === 0 ? qsTr("Covariance MC") : qsTr("Observation MC"); font.pixelSize: 18; font.weight: Font.DemiBold }
        CustomText { text: qsTr("Orbital MC"); color: AppTheme.accent; font.pixelSize: 10 }
        Item { Layout.fillWidth: true }
    }
    RowLayout {
        Layout.fillWidth: true
        Layout.topMargin: -10
        spacing: 9
        CustomIcon { name: "file-text"; color: AppTheme.muted; size: 16 }
        CustomText { text: form.monitor.method === 0 ? qsTr("Nominal + covariance") : qsTr("Observations + refit"); secondary: true; font.pixelSize: 11 }
    }
    Heading { label: form.monitor.intervalMode === 0 ? qsTr("Calculation window") : qsTr("Search interval") }
    GridLayout {
        Layout.fillWidth: true
        columns: width < 460 ? 1 : 2
        columnSpacing: 20
        rowSpacing: 18
        Field {
            label: qsTr("Event date")
            CustomComboBox {
                objectName: "risk-date-source"
                Layout.fillWidth: true
                model: [qsTr("Set date"), qsTr("Date range")]
                currentIndex: form.monitor.intervalMode
                Accessible.name: qsTr("Event date source")
                onActivated: form.monitor.intervalMode = currentIndex
            }
        }
        Field {
            visible: form.monitor.intervalMode === 0
            label: qsTr("Window centre, UTC")
            CustomTextField {
                objectName: "risk-center"
                Layout.fillWidth: true
                text: form.monitor.centerUtc
                placeholderText: "YYYY-MM-DD HH:MM:SS"
                font.family: UI.monoFont
                font.pixelSize: 12
                Accessible.name: qsTr("Window centre UTC")
                onTextEdited: form.monitor.centerUtc = text
            }
        }
        Field {
            visible: form.monitor.intervalMode === 0
            label: qsTr("Half-window, hours")
            CustomNumberField {
                objectName: "risk-half-window"
                Layout.fillWidth: true
                value: form.monitor.halfWindowHours
                minimumValue: 0
                exclusiveMinimum: true
                Accessible.name: qsTr("Window half-width in hours")
                onValueEdited: value => form.monitor.halfWindowHours = value
            }
        }
        Field {
            visible: form.monitor.intervalMode === 1
            label: qsTr("Start, UTC")
            CustomTextField {
                objectName: "risk-range-start"
                Layout.fillWidth: true
                text: form.monitor.rangeStartUtc
                placeholderText: "YYYY-MM-DD HH:MM:SS"
                font.family: UI.monoFont
                font.pixelSize: 12
                Accessible.name: qsTr("Search interval start UTC")
                onTextEdited: form.monitor.rangeStartUtc = text
            }
        }
        Field {
            visible: form.monitor.intervalMode === 1
            label: qsTr("End, UTC")
            CustomTextField {
                objectName: "risk-range-end"
                Layout.fillWidth: true
                text: form.monitor.rangeEndUtc
                placeholderText: "YYYY-MM-DD HH:MM:SS"
                font.family: UI.monoFont
                font.pixelSize: 12
                Accessible.name: qsTr("Search interval end UTC")
                onTextEdited: form.monitor.rangeEndUtc = text
            }
        }
        Field {
            visible: form.monitor.intervalMode === 1
            label: qsTr("Duration, years")
            CustomNumberField {
                objectName: "risk-range-years"
                Layout.fillWidth: true
                value: form.monitor.rangeYears
                minimumValue: 0
                exclusiveMinimum: true
                Accessible.name: qsTr("Search duration in Julian years")
                onValueEdited: value => form.monitor.rangeYears = value
            }
        }
        Field {
            label: qsTr("Target")
            CustomComboBox {
                objectName: "risk-target"
                Layout.fillWidth: true
                model: [qsTr("Earth"), qsTr("Moon"), qsTr("Mars")]
                currentIndex: form.monitor.target
                Accessible.name: qsTr("Target body")
                onActivated: form.monitor.target = currentIndex
            }
        }
        Field {
            label: qsTr("Boundary")
            CustomComboBox {
                objectName: "risk-boundary"
                Layout.fillWidth: true
                model: form.monitor.target === 0 ? [qsTr("Surface"), qsTr("Atmospheric entry \u00b7 100 km")] : [qsTr("Surface")]
                currentIndex: form.monitor.boundary
                Accessible.name: qsTr("Collision boundary")
                onActivated: form.monitor.boundary = currentIndex
            }
        }
    }
    Heading { label: qsTr("Sampling") }
    GridLayout {
        Layout.fillWidth: true
        columns: width < 460 ? 1 : 2
        columnSpacing: 20
        rowSpacing: 18
        Field {
            Layout.preferredWidth: 1
            label: form.monitor.method === 0 ? qsTr("Covariance samples") : qsTr("Observation refits")
            CustomTextField {
                objectName: "risk-samples"
                Layout.fillWidth: true
                text: form.monitor.samples
                inputMethodHints: Qt.ImhDigitsOnly
                font.family: UI.monoFont
                Accessible.name: qsTr("Sample count")
                onTextEdited: form.monitor.samples = text
            }
        }
        Field {
            Layout.preferredWidth: 1
            label: qsTr("Seed")
            CustomTextField {
                objectName: "risk-seed"
                Layout.fillWidth: true
                text: form.monitor.seed
                inputMethodHints: Qt.ImhDigitsOnly
                font.family: UI.monoFont
                Accessible.name: qsTr("Random seed")
                onTextEdited: form.monitor.seed = text
            }
        }
    }
    KeyValueList {
        visible: form.monitor.method === 1
        Layout.fillWidth: true
        rows: [[qsTr("Noise"), qsTr("Accepted measurement uncertainties")],
               [qsTr("Observation mask"), qsTr("Fixed, accepted solution")],
               [qsTr("Refit model"), qsTr("Accepted force model")]]
    }
    Heading { label: qsTr("Output") }
    RowLayout {
        Layout.fillWidth: true
        spacing: 12
        CustomIcon { name: "file-text"; color: AppTheme.blue; size: 18 }
        ColumnLayout {
            Layout.fillWidth: true
            spacing: 7
            CustomText { text: qsTr("Numerical archive + HTML report"); font.pixelSize: 11; font.weight: Font.DemiBold; Layout.fillWidth: true; wrapMode: Text.Wrap }
            CustomText { text: "JSON / OEF / RWO"; font.pixelSize: 10; secondary: true }
        }
    }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line; Layout.topMargin: 6; Layout.bottomMargin: 6 }
    CustomText {
        objectName: "risk-form-error"
        Layout.fillWidth: true
        visible: text.length > 0
        text: form.monitor.formError
        color: form.monitor.intervalMode === 1 || form.monitor.centerUtc.length ? AppTheme.amber : AppTheme.muted
        wrapMode: Text.Wrap
        elide: Text.ElideNone
        font.pixelSize: 11
    }
}
