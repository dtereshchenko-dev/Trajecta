import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI

ScrollBar {
    id: control
    InteractionCursor { enabled: control.enabled && control.interactive && control.size < 1 }
}
