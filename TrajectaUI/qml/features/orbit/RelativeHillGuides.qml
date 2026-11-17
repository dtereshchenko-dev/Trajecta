import QtQuick
import QtQuick3D
import TrajectaUI

Node {
    id: guides
    required property RelativeMotionModel model
    required property PerspectiveCamera camera
    property bool showHill: true

    Node {
        objectName: "relative-hill-billboard"
        rotation: guides.camera.sceneRotation
        Model {
            visible: guides.showHill && guides.model.relativeAvailable
            geometry: PolylineGeometry {
                objectName: "relative-hill-dashes"
                independentSegments: true
                points: {
                    const radius = guides.model.relativeHillRadius
                    if (radius <= 0) return []
                    const result = []
                    const count = 128
                    for (let i = 0; i < count; ++i) {
                        const start = 2 * Math.PI * i / count
                        const end = 2 * Math.PI * (i + 0.6) / count
                        result.push(Qt.vector3d(radius * Math.cos(start), radius * Math.sin(start), 0))
                        result.push(Qt.vector3d(radius * Math.cos(end), radius * Math.sin(end), 0))
                    }
                    return result
                }
            }
            materials: PrincipledMaterial { lighting: PrincipledMaterial.NoLighting; baseColor: AppTheme.muted }
        }
    }
    OrbitBody { objectName: "relative-l1"; visible: guides.showHill && guides.model.relativeAvailable; position: guides.model.relativeSunDirection.times(guides.model.relativeHillRadius); radius: guides.camera.z * 0.0015; color: AppTheme.amber }
    OrbitBody { objectName: "relative-l2"; visible: guides.showHill && guides.model.relativeAvailable; position: guides.model.relativeSunDirection.times(-guides.model.relativeHillRadius); radius: guides.camera.z * 0.0015; color: AppTheme.amber }
    Model {
        visible: guides.showHill && guides.model.relativeAvailable
        geometry: PolylineGeometry {
            objectName: "relative-lagrange-guides"
            independentSegments: true
            points: {
                const direction = guides.model.relativeSunDirection
                const radius = guides.model.relativeHillRadius
                if (radius <= 0) return []
                const points = []
                const count = 24
                for (const sign of [-1, 1]) {
                    for (let i = 0; i < count; ++i) {
                        points.push(direction.times(sign * radius * i / count))
                        points.push(direction.times(sign * radius * (i + 0.55) / count))
                    }
                }
                return points
            }
        }
        materials: PrincipledMaterial { lighting: PrincipledMaterial.NoLighting; baseColor: AppTheme.muted }
    }
    Model {
        visible: guides.showHill && guides.model.relativeAvailable
        geometry: PolylineGeometry {
            objectName: "relative-lagrange-rings"
            independentSegments: true
            points: {
                const radius = guides.model.relativeHillRadius
                const direction = guides.model.relativeSunDirection.normalized()
                if (radius <= 0 || direction.length() === 0) return []
                const horizontal = Math.sqrt(direction.x * direction.x + direction.y * direction.y)
                const u = horizontal > 1e-6
                          ? Qt.vector3d(-direction.y / horizontal, direction.x / horizontal, 0)
                          : Qt.vector3d(1, 0, 0)
                const v = Qt.vector3d(direction.y * u.z - direction.z * u.y,
                                      direction.z * u.x - direction.x * u.z,
                                      direction.x * u.y - direction.y * u.x)
                const ringRadius = radius * 0.1
                const points = []
                const count = 32
                for (const sign of [-1, 1]) {
                    const center = direction.times(sign * radius)
                    const at = angle => center.plus(u.times(ringRadius * Math.cos(angle)))
                                              .plus(v.times(ringRadius * Math.sin(angle)))
                    for (let i = 0; i < count; ++i) {
                        const start = 2 * Math.PI * i / count
                        const step = 2 * Math.PI / count
                        points.push(at(start), at(start + step * 0.5))
                        points.push(at(start + step * 0.7), at(start + step * 0.8))
                    }
                }
                return points
            }
        }
        materials: PrincipledMaterial { lighting: PrincipledMaterial.NoLighting; baseColor: AppTheme.amber }
    }
}
