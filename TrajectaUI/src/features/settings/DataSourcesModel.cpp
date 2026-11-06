#include "DataSourcesModel.h"
#include "settings/FitSettingsAggregator.h"
#include "orbital_core/Types.h"

#include <QDir>
#include <QDateTime>
#include <QFileInfo>
#include <QVariantMap>
#include <QtConcurrent/QtConcurrentRun>

namespace {
	QString title(const QString& id)
	{
		if (id == "planetary")
			return QObject::tr("Planetary ephemeris");
		if (id == "asteroids")
			return QObject::tr("Massive asteroids");
		if (id == "finals")
			return QObject::tr("Earth orientation: finals2000A");
		if (id == "c04")
			return QObject::tr("Earth orientation: C04");
		if (id == "delta_t")
			return QObject::tr("Historical time: TT-UT1");
		if (id == "weights")
			return QObject::tr("Astrometric weights");
		if (id == "debias_fcct14")
			return QObject::tr("Astrometric debias: FCCT14");
		if (id == "debias_efcc18")
			return QObject::tr("Astrometric debias: EFCC18");
		if (id == "observatories")
			return QObject::tr("Optical observatories");
		if (id == "radar")
			return QObject::tr("Radar observatories");
		return id;
	}

	QString impact(const QString& id)
	{
		if (id == "planetary")
			return QObject::tr("Needed for planetary ephemerides and N-body calculations.");
		if (id == "asteroids")
			return QObject::tr("Needed when massive-asteroid perturbations are enabled.");
		if (id == "finals")
			return QObject::tr("Used for UT1; other time data or an approximation may be used outside its coverage.");
		if (id == "c04")
			return QObject::tr("Used for the Earth pole in the gravity model; missing data uses a nominal pole.");
		if (id == "delta_t")
			return QObject::tr("Used for historical TT–UT1 conversion when finals data do not cover the epoch.");
		if (id == "weights")
			return QObject::tr("Needed for VFCC17 weights; missing rules fall back to 1 arcsec.");
		if (id == "debias_fcct14" || id == "debias_efcc18")
			return QObject::tr("Needed only when the corresponding catalogue debias is enabled.");
		if (id == "observatories")
			return QObject::tr("Needed to load observations and resolve observatory codes.");
		if (id == "radar")
			return QObject::tr("Needed for radar observations; optical observations do not require it.");
		return {};
	}
}  // namespace

DataSourcesModel::DataSourcesModel(FitSettingsAggregator* settings, QObject* parent)
    : QObject(parent), m_settings(settings)
{
	if (m_settings)
		connect(m_settings, &FitSettingsAggregator::anyChanged, this, [this] {
			if (m_checked) {
				updateReadiness();
				emit changed();
			}
		});
	connect(&m_sb441Download, &Sb441Download::installed, this, [this] {
		m_restartRequired = true;
		refresh();
	});
	connect(&m_observatoryUpdate, &ObservatoryUpdate::installed, this, [this] {
		m_restartRequired = true;
		refresh();
	});
	connect(&m_observatoryUpdate, &ObservatoryUpdate::changed, this, [this] {
		if (m_checked) updateReadiness();
		emit changed();
	});
	m_debounce.setSingleShot(true);
	m_debounce.setInterval(750);
	connect(&m_debounce, &QTimer::timeout, this, &DataSourcesModel::refresh);
	connect(&m_files, &QFileSystemWatcher::fileChanged, this, [this] {
		if (m_checked) {
			m_restartRequired = true;
			emit changed();
		}
		m_debounce.start();
	});
	connect(&m_files, &QFileSystemWatcher::directoryChanged, this, [this] { m_debounce.start(); });
	connect(&m_watcher, &QFutureWatcherBase::finished, this, [this] {
		const auto previousEntries = m_entries;
		m_entries.clear();
		m_observationsAvailable = false;
		for (const auto& source : m_watcher.result()) {
			const auto id = QString::fromStdString(source.id);
			const auto path = QString::fromStdWString(source.path.wstring());
			const QFileInfo file(path);
			const bool onDisk = file.isFile() && file.isReadable();
			if (id == "observatories")
				m_observationsAvailable = source.active && source.valid;
			QStringList details{file.fileName()};
			if (source.records)
				details << tr("%1 records").arg(source.records);
			if (source.catalogs)
				details << tr("%1 catalogues").arg(source.catalogs);
			if (source.firstMjdUtc && source.lastMjdUtc)
				details << tr("MJD %1–%2 UTC").arg(source.firstMjdUtc, 0, 'f', 0).arg(source.lastMjdUtc, 0, 'f', 0);
			if (file.isFile())
				details << (file.size() < 1048576 ? tr("%1 KiB").arg(file.size() / 1024.0, 0, 'f', 1)
				                                  : tr("%1 MiB").arg(file.size() / 1048576.0, 0, 'f', 1));
			if (!source.error.empty())
				details << QString::fromStdString(source.error);
			m_entries.append(
			    QVariantMap{{"id", id},
			                {"title", title(id)},
			                {"path", path},
			                {"detail", details.join(QStringLiteral("  ·  "))},
			                {"valid", source.valid},
			                {"impact", impact(id)},
			                {"active", source.active},
			                {"onDisk", onDisk},
			                {"sizeBytes", file.size()},
			                {"modifiedMs", file.lastModified().toMSecsSinceEpoch()},
			                {"status", !onDisk         ? (source.active ? tr("Loaded; file missing") : tr("Missing"))
			                           : source.active ? tr("Connected")
			                           : source.valid  ? tr("Verified on disk")
			                                           : tr("Invalid")}});
		}
		if (m_checked) {
			for (const auto& value : m_entries) {
				const auto entry = value.toMap();
				for (const auto& oldValue : previousEntries) {
					const auto old = oldValue.toMap();
					if (old.value("id") != entry.value("id"))
						continue;
					if (old.value("onDisk") != entry.value("onDisk") || old.value("valid") != entry.value("valid") ||
					    old.value("sizeBytes") != entry.value("sizeBytes") ||
					    old.value("modifiedMs") != entry.value("modifiedMs"))
						m_restartRequired = true;
					break;
				}
			}
		}
		m_checked = true;
		m_running = false;
		updateReadiness();
		updateWatchPaths();
		emit changed();
		if (m_refreshPending) {
			m_refreshPending = false;
			refresh();
		}
	});
	updateWatchPaths();
}

void DataSourcesModel::updateReadiness()
{
	QStringList unavailable, selected;
	const auto* pert = m_settings ? m_settings->perturbations() : nullptr;
	const auto* weight = m_settings ? m_settings->weighting() : nullptr;
	for (const auto& value : m_entries) {
		const auto entry = value.toMap();
		const auto id = entry.value("id").toString();
		const auto name = entry.value("title").toString();
		const bool configured =
		    (id == "planetary" && pert && pert->enabled()) ||
		    (id == "asteroids" && pert && pert->enabled() && pert->useSb441()) ||
		    (id == "weights" && weight && weight->useWeights() && weight->weightProfile() == 1) ||
		    (id == "debias_efcc18" && weight && weight->debiasEnabled());
		const bool availableOnDisk = entry.value("valid").toBool() && entry.value("onDisk").toBool();
		if (!availableOnDisk)
			unavailable.append(name);
		if (configured && !entry.value("active").toBool() && (!availableOnDisk || id == "asteroids"))
			selected.append(name);
	}
	if (!m_observationsAvailable)
		selected.prepend(tr("Optical observatories"));
	m_unavailableCount = unavailable.size();
	m_selectedUnavailableCount = selected.size();
	const auto format = [this](const QStringList& names, const QString& one, const QString& many) {
		return names.size() == 1 ? one.arg(names.first()) : many.arg(names.first()).arg(names.size() - 1);
	};
	if (!selected.isEmpty())
		m_readinessMessage =
		    format(selected, tr("Selected data unavailable: %1"), tr("Selected data unavailable: %1 and %2 more"));
	else if (!unavailable.isEmpty())
		m_readinessMessage =
		    format(unavailable, tr("Ready with limited data: %1"), tr("Ready with limited data: %1 and %2 more"));
	else if (m_restartRequired)
		m_readinessMessage = tr("Data changed; restart required");
	else if (m_observatoryUpdate.availableCount() > 0)
		m_readinessMessage = tr("Observatory updates available");
	else
		m_readinessMessage = tr("Ready");
}

DataSourcesModel::~DataSourcesModel() { m_watcher.waitForFinished(); }

void DataSourcesModel::refresh()
{
	if (checking()) {
		m_refreshPending = true;
		return;
	}
	m_running = true;
	const auto installed = std::filesystem::u8path(Sb441Download::installedPath().toUtf8().toStdString());
	m_watcher.setFuture(QtConcurrent::run([installed] { return orbital_core::inspectRuntimeData(installed); }));
	emit changed();
}

void DataSourcesModel::updateWatchPaths()
{
	const QDir base(QString::fromStdWString(getExecutableDir().wstring()));
	const auto downloaded = QFileInfo(Sb441Download::installedPath()).absolutePath();
	const auto observatories = QString::fromStdWString(orbital_core::observatoryDataDirectory().wstring());
	const QStringList directories{base.absolutePath(),          base.filePath("earth_orientation"),
	                              base.filePath("weights"),     base.filePath("debias"),
	                              base.filePath("debias_fcct14"), base.filePath("debias_2018"), downloaded,
	                              observatories};
	QStringList wantedDirectories;
	for (const auto& directory : directories)
		if (QFileInfo(directory).isDir())
			wantedDirectories.append(directory);
	for (const auto& directory : m_files.directories())
		if (!wantedDirectories.contains(directory))
			m_files.removePath(directory);
	for (const auto& directory : wantedDirectories)
		if (!m_files.directories().contains(directory))
			m_files.addPath(directory);
	QStringList wantedFiles;
	for (const auto& value : m_entries) {
		const auto path = value.toMap().value("path").toString();
		if (QFileInfo(path).isFile())
			wantedFiles.append(path);
	}
	const auto corrections = base.filePath("earth_orientation/eopc04.cnc");
	if (QFileInfo(corrections).isFile())
		wantedFiles.append(corrections);
	for (const auto& path : m_files.files())
		if (!wantedFiles.contains(path))
			m_files.removePath(path);
	for (const auto& path : wantedFiles)
		if (!m_files.files().contains(path))
			m_files.addPath(path);
}
