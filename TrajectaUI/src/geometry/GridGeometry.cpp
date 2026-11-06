#include "geometry/GridGeometry.h"
#include <QVector3D>

GridGeometry::GridGeometry(QQuick3DObject* parent) : QQuick3DGeometry(parent) { rebuild(); }

void GridGeometry::rebuild()
{
	clear();

	const float half = m_size / 2.0f;
	const float step = m_size / m_divisions;
	const int lines = m_divisions + 1;

	const int totalVerts = lines * 2 * 2;

	QByteArray vertexData;
	vertexData.resize(totalVerts * 3 * sizeof(float));
	auto* vptr = reinterpret_cast<float*>(vertexData.data());

	QByteArray indexData;
	indexData.resize(lines * 2 * 2 * sizeof(quint32));
	auto* iptr = reinterpret_cast<quint32*>(indexData.data());

	quint32 idx = 0;

	for (int i = 0; i < lines; ++i) {
		float x = -half + i * step;
		*vptr++ = x;
		*vptr++ = 0;
		*vptr++ = -half;
		*vptr++ = x;
		*vptr++ = 0;
		*vptr++ = half;
		*iptr++ = idx;
		*iptr++ = idx + 1;
		idx += 2;
	}

	for (int i = 0; i < lines; ++i) {
		float z = -half + i * step;
		*vptr++ = -half;
		*vptr++ = 0;
		*vptr++ = z;
		*vptr++ = half;
		*vptr++ = 0;
		*vptr++ = z;
		*iptr++ = idx;
		*iptr++ = idx + 1;
		idx += 2;
	}

	setVertexData(vertexData);
	setIndexData(indexData);
	setStride(3 * sizeof(float));

	addAttribute(QQuick3DGeometry::Attribute::PositionSemantic, 0, QQuick3DGeometry::Attribute::F32Type);
	addAttribute(QQuick3DGeometry::Attribute::IndexSemantic, 0, QQuick3DGeometry::Attribute::U32Type);

	setPrimitiveType(QQuick3DGeometry::PrimitiveType::Lines);
	setBounds(QVector3D(-half, 0, -half), QVector3D(half, 0, half));

	update();
}
