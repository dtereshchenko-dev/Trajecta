import QtQuick
import QtQuick.Layouts
import QtQuick3D
import TrajectaUI

RowLayout {
    id: indicator
    objectName: "relative-sun-indicator"
    required property RelativeMotionModel model
    required property PerspectiveCamera camera
    spacing: 6
    visible: model.relativeAvailable

    View3D {
        Layout.preferredWidth: 28
        Layout.preferredHeight: 28
        camera: indicatorCamera
        environment: SceneEnvironment { backgroundMode: SceneEnvironment.Transparent; antialiasingMode: SceneEnvironment.MSAA }
        Node {
            objectName: "relative-sun-camera-origin"
            rotation: indicator.camera.sceneRotation
            OrthographicCamera { id: indicatorCamera; z: 100; clipNear: 1; clipFar: 200 }
        }
        Node {
            objectName: "relative-sun-arrow"
            readonly property vector3d direction: indicator.model.relativeSunDirection.normalized()
            rotation: {
                const d = direction
                if (d.y < -0.999999) return Qt.quaternion(0, 1, 0, 0)
                const length = Math.sqrt(2 * (1 + d.y))
                return Qt.quaternion(length / 2, d.z / length, 0, -d.x / length)
            }
            Model {
                source: "#Cylinder"
                y: -3
                scale: Qt.vector3d(0.018, 0.14, 0.018)
                materials: PrincipledMaterial { lighting: PrincipledMaterial.NoLighting; baseColor: AppTheme.amber }
            }
            Model {
                source: "#Cone"
                y: 7
                scale: Qt.vector3d(0.055, 0.07, 0.055)
                materials: PrincipledMaterial { lighting: PrincipledMaterial.NoLighting; baseColor: AppTheme.amber }
            }
        }
    }
    CustomText { text: qsTr("Sun"); color: AppTheme.amber; font.pixelSize: 11 }
}
