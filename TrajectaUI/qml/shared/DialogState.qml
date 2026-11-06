pragma Singleton
import QtQuick

QtObject {
    property var dialogs: []

    function update(dialog, window, active) {
        const next = dialogs.filter(entry => entry.dialog !== dialog)
        if (active && window) next.push({ dialog: dialog, window: window })
        dialogs = next
    }

    function activeFor(window) {
        return dialogs.some(entry => entry.window === window)
    }
}
