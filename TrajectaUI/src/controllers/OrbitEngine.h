#pragma once

#include <QObject>
#include <QThread>
#include "controllers/OrbitWorker.h"
#include "settings/FitSettingsAggregator.h"

class OrbitEngine : public QObject {
	Q_OBJECT

   public:
	explicit OrbitEngine(FitSettingsAggregator& settings, QObject* parent = nullptr);
	~OrbitEngine() override;

	QString errorMessage() const { return m_errorMessage; }

	QString progressStatus() const { return m_progressStatus; }

	void loadFile(const QString& path);
	void prepareFile(const QString& path);
	void recalculate(const orbital_core::ObservationSet& observations,
	                 std::vector<orbital_core::ObservationOverride> overrides);
	void downloadGaiaFpr(orbital_core::ObservationSet observations,
	                     std::vector<orbital_core::ObservationOverride> overrides,
	                     std::optional<double> photocenterRadiusKm);
	void cancel();
	OptimizationSettings currentSettings() const;

   signals:
	void errorMessageChanged();
	void progressStatusChanged();
	void solutionReady(const orbital_core::FitPipelineResult& result, const QString& directory);
	void observationsReady(std::shared_ptr<orbital_core::ObservationSet> observations);

   private:
	void submit(OrbitFitInput input);

	QThread m_workerThread;
	OrbitWorker* m_worker{};
	std::stop_source m_fitCancel;
	bool m_fitRunning{};
	QString m_errorMessage;
	QString m_progressStatus;
	FitSettingsAggregator& m_settings;
};
