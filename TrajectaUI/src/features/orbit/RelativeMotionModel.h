#pragma once

#include <QQmlEngine>
#include <QVector3D>
#include <orbital_core/RelativeMotion.h>
#include <optional>
#include <thread>

namespace orbital_core {
	struct FitPipelineResult;
}

class RelativeMotionModel : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Owned by the orbit scene")
	Q_PROPERTY(bool relativeBusy READ relativeBusy NOTIFY statusChanged)
	Q_PROPERTY(bool relativeAvailable READ relativeAvailable NOTIFY relativeChanged)
	Q_PROPERTY(bool rotatingFrame READ rotatingFrame WRITE setRotatingFrame NOTIFY relativeChanged)
	Q_PROPERTY(QVariantList relativePath READ relativePath NOTIFY relativeChanged)
	Q_PROPERTY(QVariantList relativeBoundSegments READ relativeBoundSegments NOTIFY relativeChanged)
	Q_PROPERTY(QVariantList relativeUnboundSegments READ relativeUnboundSegments NOTIFY relativeChanged)
	Q_PROPERTY(QVariantList relativeDistantNegativeSegments READ relativeDistantNegativeSegments NOTIFY relativeChanged)
	Q_PROPERTY(QVariantList relativeMoonPath READ relativeMoonPath NOTIFY relativeChanged)
	Q_PROPERTY(QVariantList relativeCaptures READ relativeCaptures NOTIFY relativeChanged)
	Q_PROPERTY(QVariantList relativeCaptureMarkers READ relativeCaptureMarkers NOTIFY relativeChanged)
	Q_PROPERTY(QVector3D relativeMoon READ relativeMoon NOTIFY positionChanged)
	Q_PROPERTY(QVector3D relativeSunDirection READ relativeSunDirection NOTIFY positionChanged)
	Q_PROPERTY(double relativeHillRadius READ relativeHillRadius NOTIFY positionChanged)
	Q_PROPERTY(double relativePlanetRadius READ relativePlanetRadius NOTIFY relativeChanged)
	Q_PROPERTY(double relativeDistanceKm READ relativeDistanceKm NOTIFY positionChanged)
	Q_PROPERTY(double relativeEnergy READ relativeEnergy NOTIFY positionChanged)
	Q_PROPERTY(QString relativeError READ relativeError NOTIFY statusChanged)
	Q_PROPERTY(QString relativeCenter READ relativeCenter NOTIFY relativeChanged)
	Q_PROPERTY(QString relativeStartDate READ relativeStartDate WRITE setRelativeStartDate NOTIFY settingsChanged)
	Q_PROPERTY(QString relativeEndDate READ relativeEndDate WRITE setRelativeEndDate NOTIFY settingsChanged)
	Q_PROPERTY(double relativeStepDays READ relativeStepDays WRITE setRelativeStepDays NOTIFY settingsChanged)
	Q_PROPERTY(int relativePlanetIndex READ relativePlanetIndex WRITE setRelativePlanetIndex NOTIFY settingsChanged)
	Q_PROPERTY(QString relativeTimeScale READ relativeTimeScale WRITE setRelativeTimeScale NOTIFY settingsChanged)
	Q_PROPERTY(QString relativeResultTimeScale READ relativeResultTimeScale NOTIFY relativeChanged)

   public:
	explicit RelativeMotionModel(QObject* parent = nullptr);
	~RelativeMotionModel() override;
	void reset();
	void setSolution(const orbital_core::FitPipelineResult& result);
	void setPosition(double position);

	QVector3D objectPosition() const { return m_objectPosition; }

	double firstMjdTdb() const { return m_relativeSamples.front().mjdTdb; }

	double spanDays() const { return m_relativeSamples.back().mjdTdb - m_relativeSamples.front().mjdTdb; }

	double frameRadius() const { return m_relativeRadius; }

	std::optional<double> captureStartEpoch(int index) const;

	bool relativeBusy() const { return m_relativeBusy; }

	bool relativeAvailable() const { return m_relativeSamples.size() >= 2; }

	bool rotatingFrame() const { return m_rotatingFrame; }

	void setRotatingFrame(bool value);

	QVariantList relativePath() const { return m_relativePath; }

	QVariantList relativeBoundSegments() const { return m_relativeBoundSegments; }

	QVariantList relativeUnboundSegments() const { return m_relativeUnboundSegments; }

	QVariantList relativeDistantNegativeSegments() const { return m_relativeDistantNegativeSegments; }

	QVariantList relativeMoonPath() const { return m_relativeMoonPath; }

	QVariantList relativeCaptures() const { return m_relativeCaptures; }

	QVariantList relativeCaptureMarkers() const { return m_relativeCaptureMarkers; }

	QVector3D relativeMoon() const { return m_relativeMoon; }

	QVector3D relativeSunDirection() const;
	double relativeHillRadius() const;

	double relativePlanetRadius() const { return m_relativePlanetRadius; }

	double relativeDistanceKm() const;
	double relativeEnergy() const;

	QString relativeError() const { return m_relativeError; }

	QString relativeCenter() const { return m_relativeCenter; }

	QString relativeStartDate() const { return m_relativeStartDate; }

	void setRelativeStartDate(const QString& value);

	QString relativeEndDate() const { return m_relativeEndDate; }

	void setRelativeEndDate(const QString& value);

	double relativeStepDays() const { return m_relativeStepDays; }

	void setRelativeStepDays(double value);

	int relativePlanetIndex() const { return m_relativePlanetIndex; }

	void setRelativePlanetIndex(int value);

	QString relativeTimeScale() const { return m_relativeTimeScale; }

	void setRelativeTimeScale(const QString& value);

	QString relativeResultTimeScale() const { return m_relativeResultTimeScale; }

	Q_INVOKABLE void calculateRelative(const QString& startUtc, const QString& endUtc, int planetIndex,
	                                   double stepDays);
	Q_INVOKABLE void cancelRelative();

   signals:
	void relativeChanged();
	void settingsChanged();
	void statusChanged();
	void positionChanged();
	void resultReady(bool accepted);

   private:
	void updateRelativePosition();
	void rebuildRelativePath();
	const orbital_core::RelativeMotionSample* currentRelativeSample() const;
	std::optional<orbital_core::EphemerisSolution> m_relativeSolution;
	std::vector<orbital_core::RelativeMotionSample> m_relativeSamples;
	std::vector<orbital_core::TemporaryCaptureInterval> m_captureIntervals;
	QVariantList m_relativePath, m_relativeMoonPath;
	QVariantList m_relativeBoundSegments, m_relativeUnboundSegments;
	QVariantList m_relativeDistantNegativeSegments, m_relativeCaptures, m_relativeCaptureMarkers;
	QVector3D m_objectPosition, m_relativeMoon;
	std::jthread m_relativeWorker;
	quint64 m_relativeGeneration{};
	double m_position{}, m_relativeRadius{0.01}, m_relativePlanetRadius{}, m_relativeStepDays{0.1};
	int m_relativePlanetIndex{2};
	bool m_relativeBusy{}, m_rotatingFrame{true};
	QString m_relativeError, m_relativeCenter, m_relativeStartDate, m_relativeEndDate;
	QString m_relativeTimeScale{"UTC"}, m_relativeResultTimeScale{"UTC"};
};
