import QtQuick
import QtQuick3D
import TrajectaUI

Model {
    id: path
    property real a: 1
    property real e: 0
    property real q: 0
    property real inclination: 0
    property real node: 0
    property real peri: 0
    property color color: AppTheme.accent
    readonly property vector3d perihelionPosition: curve.perihelionPos
    geometry: OrbitGeometry {
        id: curve
        semiMajorAxis: path.a
        eccentricity: path.e
        perigeeAu: path.q
        inclination: path.inclination
        ascendingNode: path.node
        argPericentre: path.peri
        segments: 512
    }
    materials: PrincipledMaterial { lighting: PrincipledMaterial.NoLighting; baseColor: path.color }
    castsShadows: false
    receivesShadows: false
}
