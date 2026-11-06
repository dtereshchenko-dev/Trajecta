import QtQuick
import Qt5Compat.GraphicalEffects

Item {
    property bool active: false
    layer.enabled: active && GraphicsInfo.api !== GraphicsInfo.Software
    layer.effect: GaussianBlur {
        radius: Math.ceil(6 * Screen.devicePixelRatio)
        samples: radius * 2 + 1
        deviation: 2 * Screen.devicePixelRatio
        transparentBorder: false
    }
}
