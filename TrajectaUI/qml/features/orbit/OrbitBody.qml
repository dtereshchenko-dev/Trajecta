import QtQuick
import QtQuick3D
import TrajectaUI

Model {
    id: body
    property color color: AppTheme.accent
    property real radius: 0.015
    source: "#Sphere"
    scale: Qt.vector3d(radius / 50, radius / 50, radius / 50)
    materials: PrincipledMaterial {
        lighting: PrincipledMaterial.NoLighting
        baseColor: body.color
    }
    castsShadows: false
    receivesShadows: false
}
