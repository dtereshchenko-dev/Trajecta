#pragma once

#include <QObject>
#include <QQmlEngine>
#include <QFutureWatcher>
#include <QFileSystemWatcher>
#include <QTimer>
#include <QVariantList>
#include "orbital_core/RuntimeData.h"
#include "features/settings/Sb441Download.h"
#include "features/settings/ObservatoryUpdate.h"

class FitSettingsAggregator;

class DataSourcesModel : public QObject {
	Q_OBJECT
	QML_ELEMENT
	QML_UNCREATABLE("Session-owned")
	Q_PROPERTY(QVariantList entries READ entries NOTIFY changed)
	Q_PROPERTY(bool checking READ checking NOTIFY changed)
	Q_PROPERTY(bool checked READ checked NOTIFY changed)
	Q_PROPERTY(int unavailableCount READ unavailableCount NOTIFY changed)
	Q_PROPERTY(bool restartRequired READ restartRequired NOTIFY changed)
	Q_PROPERTY(bool observationsAvailable READ observationsAvailable NOTIFY changed)
	Q_PROPERTY(int selectedUnavailableCount READ selectedUnavailableCount NOTIFY changed)
	Q_PROPERTY(bool hasWarning READ hasWarning NOTIFY changed)
	Q_PROPERTY(QString readinessMessage READ readinessMessage NOTIFY changed)
	Q_PROPERTY(bool planetaryAvailable READ planetaryAvailable CONSTANT)
	Q_PROPERTY(bool asteroidsAvailable READ asteroidsAvailable CONSTANT)
	Q_PROPERTY(Sb441Download* sb441Download READ sb441Download CONSTANT)
	Q_PROPERTY(ObservatoryUpdate* observatoryUpdate READ observatoryUpdate CONSTANT)

   public:
	explicit DataSourcesModel(FitSettingsAggregator* settings = nullptr, QObject* parent = nullptr);
	~DataSourcesModel() override;

	QVariantList entries() const { return m_entries; }

	bool checking() const { return m_running; }

	bool checked() const { return m_checked; }

	int unavailableCount() const { return m_unavailableCount; }

	bool restartRequired() const { return m_restartRequired; }

	bool observationsAvailable() const { return m_observationsAvailable; }

	int selectedUnavailableCount() const { return m_selectedUnavailableCount; }

	bool hasWarning() const
	{
		return (m_checked && (m_unavailableCount > 0 || m_selectedUnavailableCount > 0 || m_restartRequired)) ||
		       m_observatoryUpdate.availableCount() > 0;
	}

	QString readinessMessage() const { return m_readinessMessage; }

	bool planetaryAvailable() const { return orbital_core::planetaryEphemerisAvailable(); }

	bool asteroidsAvailable() const { return orbital_core::asteroidEphemerisAvailable(); }

	Sb441Download* sb441Download() { return &m_sb441Download; }

	ObservatoryUpdate* observatoryUpdate() { return &m_observatoryUpdate; }

	Q_INVOKABLE void refresh();

   signals:
	void changed();

   private:
	void updateWatchPaths();
	void updateReadiness();
	FitSettingsAggregator* m_settings;
	Sb441Download m_sb441Download;
	ObservatoryUpdate m_observatoryUpdate;
	QFutureWatcher<std::vector<orbital_core::RuntimeDataEntry>> m_watcher;
	QFileSystemWatcher m_files;
	QTimer m_debounce;
	QVariantList m_entries;
	QString m_readinessMessage;
	int m_unavailableCount = 0;
	int m_selectedUnavailableCount = 0;
	bool m_checked = false;
	bool m_refreshPending = false;
	bool m_running = false;
	bool m_restartRequired = false;
	bool m_observationsAvailable = false;
};
