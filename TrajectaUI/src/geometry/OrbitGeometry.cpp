#include "geometry/OrbitGeometry.h"
#include <algorithm>
#include <limits>
#include <cmath>
#include "orbital_core/Types.h"
#include <QByteArray>

OrbitGeometry::OrbitGeometry(QQuick3DObject* parent) : QQuick3DGeometry(parent) { rebuild(); }

QVector3D OrbitGeometry::orbitalToEcliptic(double x_orb, double y_orb) const
{
	const double w = Convertors::degreesToRadians(m_omega);
	const double O = Convertors::degreesToRadians(m_Omega);
	const double i = Convertors::degreesToRadians(m_i);

	const double cw = std::cos(w), sw = std::sin(w);
	const double cO = std::cos(O), sO = std::sin(O);
	const double ci = std::cos(i), si = std::sin(i);

	double x = (cw * cO - sw * sO * ci) * x_orb + (-sw * cO - cw * sO * ci) * y_orb;
	double y = (cw * sO + sw * cO * ci) * x_orb + (-sw * sO + cw * cO * ci) * y_orb;
	double z = (sw * si) * x_orb + (cw * si) * y_orb;

	return QVector3D(static_cast<float>(x), static_cast<float>(z), static_cast<float>(-y));
}

void OrbitGeometry::rebuild()
{
	clear();
	m_perihelionPos = QVector3D();
	setBounds(QVector3D(), QVector3D());

	const double e = m_e;
	const double a = m_a;
	if (!std::isfinite(e) || e < 0.0 || !std::isfinite(m_q) || m_q < 0.0 || !std::isfinite(m_i) ||
	    !std::isfinite(m_Omega) || !std::isfinite(m_omega) || m_segments < 2 || (m_q == 0.0 && !std::isfinite(a))) {
		update();
		return;
	}

	const double p = m_q > 0.0 ? m_q * (1.0 + e) : std::abs(a) * std::abs(e * e - 1.0);
	if (!std::isfinite(p) || p <= 0.0) {
		update();
		return;
	}

	bool isElliptical = (e < 1.0);

	double nuMin, nuMax;
	if (isElliptical) {
		nuMin = 0.0;
		nuMax = 2.0 * Constants::PI;
	}

	else {
		double nuLim = std::acos(-1.0 / e) * 0.95;
		nuMin = -nuLim;
		nuMax = nuLim;
	}

	const int N = m_segments;
	const int numPoints = isElliptical ? N + 1 : N;
	const double step = (nuMax - nuMin) / N;

	QByteArray vertexData;
	vertexData.resize(numPoints * 3 * sizeof(float));
	auto* vptr = reinterpret_cast<float*>(vertexData.data());

	float minDist = std::numeric_limits<float>::max();
	QVector3D boundsMin(minDist, minDist, minDist);
	QVector3D boundsMax(-minDist, -minDist, -minDist);

	for (int i = 0; i < numPoints; ++i) {
		double nu = nuMin + i * step;
		double r = p / (1.0 + e * std::cos(nu));
		QVector3D pt = orbitalToEcliptic(r * std::cos(nu), r * std::sin(nu));
		if (!std::isfinite(pt.x()) || !std::isfinite(pt.y()) || !std::isfinite(pt.z())) {
			m_perihelionPos = QVector3D();
			update();
			return;
		}
		boundsMin.setX(std::min(boundsMin.x(), pt.x()));
		boundsMin.setY(std::min(boundsMin.y(), pt.y()));
		boundsMin.setZ(std::min(boundsMin.z(), pt.z()));
		boundsMax.setX(std::max(boundsMax.x(), pt.x()));
		boundsMax.setY(std::max(boundsMax.y(), pt.y()));
		boundsMax.setZ(std::max(boundsMax.z(), pt.z()));
		*vptr++ = pt.x();
		*vptr++ = pt.y();
		*vptr++ = pt.z();

		float dist = pt.length();
		if (dist < minDist) {
			minDist = dist;
			m_perihelionPos = pt;
		}
	}

	int numSegments = isElliptical ? N : N - 1;
	QByteArray indexData;
	indexData.resize(numSegments * 2 * sizeof(quint32));
	auto* iptr = reinterpret_cast<quint32*>(indexData.data());

	for (int i = 0; i < numSegments; ++i) {
		*iptr++ = static_cast<quint32>(i);
		*iptr++ = static_cast<quint32>((i + 1) % numPoints);
	}

	setVertexData(vertexData);
	setIndexData(indexData);
	setStride(3 * sizeof(float));

	addAttribute(QQuick3DGeometry::Attribute::PositionSemantic, 0, QQuick3DGeometry::Attribute::F32Type);
	addAttribute(QQuick3DGeometry::Attribute::IndexSemantic, 0, QQuick3DGeometry::Attribute::U32Type);

	setPrimitiveType(QQuick3DGeometry::PrimitiveType::Lines);

	setBounds(boundsMin, boundsMax);

	update();
}
