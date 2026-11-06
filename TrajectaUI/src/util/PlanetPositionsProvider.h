#pragma once

#include <QObject>
#include <QVariantList>

class PlanetPositionsProvider : public QObject {
	Q_OBJECT
	Q_PROPERTY(QVariantList positions READ positions NOTIFY positionsChanged)

   public:
	explicit PlanetPositionsProvider(QObject* parent = nullptr) : QObject(parent) {}

	void setEpoch(double jd);

	QVariantList positions() const { return m_positions; }

   signals:
	void positionsChanged();

   private:
	QVariantList m_positions;
};
