#include "PolylineGeometry.h"

#include <QVector3D>
#include <algorithm>
#include <limits>

PolylineGeometry::PolylineGeometry(QQuick3DObject* parent) : QQuick3DGeometry(parent) {}

void PolylineGeometry::setPoints(const QVariantList& points)
{
	if (m_points == points)
		return;
	m_points = points;
	rebuild();
	emit pointsChanged();
}

void PolylineGeometry::setIndependentSegments(bool value)
{
	if (m_independentSegments == value)
		return;
	m_independentSegments = value;
	rebuild();
	emit independentSegmentsChanged();
}

void PolylineGeometry::rebuild()
{
	clear();
	if (m_points.size() >= 2) {
		QByteArray vertices;
		const qsizetype verticesCount = m_independentSegments ? (m_points.size() / 2) * 2 : (m_points.size() - 1) * 2;
		vertices.resize(verticesCount * 3 * sizeof(float));
		float* out = reinterpret_cast<float*>(vertices.data());
		const float inf = std::numeric_limits<float>::infinity();
		QVector3D lo(inf, inf, inf), hi(-inf, -inf, -inf);
		const auto write = [&](const QVector3D& p) {
			*out++ = p.x();
			*out++ = p.y();
			*out++ = p.z();
			lo.setX(std::min(lo.x(), p.x()));
			lo.setY(std::min(lo.y(), p.y()));
			lo.setZ(std::min(lo.z(), p.z()));
			hi.setX(std::max(hi.x(), p.x()));
			hi.setY(std::max(hi.y(), p.y()));
			hi.setZ(std::max(hi.z(), p.z()));
		};
		for (qsizetype i = 1; i < m_points.size(); i += m_independentSegments ? 2 : 1) {
			write(m_points[i - 1].value<QVector3D>());
			write(m_points[i].value<QVector3D>());
		}
		setVertexData(vertices);
		setStride(3 * sizeof(float));
		addAttribute(Attribute::PositionSemantic, 0, Attribute::F32Type);
		setPrimitiveType(PrimitiveType::Lines);
		setBounds(lo, hi);
	}
	update();
}
