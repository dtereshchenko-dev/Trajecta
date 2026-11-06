#include "ObservatoryUpdate.h"
#include "orbital_core/RuntimeData.h"
#include "orbital_core/obs/ObservatoryUpdates.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QSaveFile>
#include <QtConcurrent/QtConcurrentRun>
#include <exception>

namespace {
	QString catalogPath(const QString& id)
	{
		const auto directory = orbital_core::observatoryDataDirectory();
		if (directory.empty()) return {};
		return QDir(QString::fromStdU16String(directory.u16string()))
		    .filePath(id == "observatories" ? "OBSCODE.dat" : "RADCODE.dat");
	}

	QByteArray readCatalog(const QString& path, bool& exists)
	{
		if (path.isEmpty()) throw std::runtime_error("Observatory data directory is not configured");
		QFile file(path);
		exists = file.exists();
		if (!exists) return {};
		if (!file.open(QIODevice::ReadOnly)) throw std::runtime_error(file.errorString().toStdString());
		if (file.size() > 4 * 1024 * 1024) throw std::runtime_error("Observatory catalog is too large");
		const auto data = file.readAll();
		if (file.error() != QFileDevice::NoError) throw std::runtime_error(file.errorString().toStdString());
		return data;
	}
}

ObservatoryUpdate::ObservatoryUpdate(QObject* parent) : ObservatoryUpdate(Fetch{}, parent) {}

ObservatoryUpdate::ObservatoryUpdate(Fetch fetch, QObject* parent) : QObject(parent), m_fetch(std::move(fetch))
{
	if (!m_fetch)
		m_fetch = [](const QString& id, std::string_view current, std::stop_token token) {
			return id == "observatories"
			    ? orbital_core::obs::downloadOpticalObservatoryUpdate(current, token)
			    : orbital_core::obs::downloadRadarObservatoryUpdate(current, token);
		};
	m_optical.id = "observatories";
	m_radar.id = "radar";
	for (auto* item : {&m_optical, &m_radar}) {
		connect(&item->watcher, &QFutureWatcherBase::finished, this, [this, item] {
			item->result = item->watcher.result();
			item->checking = false;
			item->available = item->result.error.isEmpty() && item->result.changed;
			item->message = item->result.error.isEmpty() ? QString() : tr("Update check unavailable: %1").arg(item->result.error);
			emit changed();
		});
	}
}

ObservatoryUpdate::~ObservatoryUpdate()
{
	for (auto* item : {&m_optical, &m_radar}) {
		item->cancel.request_stop();
		item->watcher.waitForFinished();
	}
}

ObservatoryUpdate::Source* ObservatoryUpdate::source(const QString& id)
{
	if (id == m_optical.id) return &m_optical;
	if (id == m_radar.id) return &m_radar;
	return nullptr;
}

int ObservatoryUpdate::availableCount() const { return int(m_optical.available) + int(m_radar.available); }
QVariantMap ObservatoryUpdate::states() const
{
	QVariantMap result;
	for (const auto* item : {&m_optical, &m_radar})
		result.insert(item->id, QVariantMap{{"available", item->available}, {"checking", item->checking},
		                                     {"message", item->message}});
	return result;
}

void ObservatoryUpdate::checkSource(Source& item)
{
	if (item.checking) return;
	const auto id = item.id;
	const auto path = catalogPath(id);
	item.cancel = std::stop_source{};
	const auto token = item.cancel.get_token();
	const auto fetch = m_fetch;
	item.checking = true;
	item.available = false;
	item.message.clear();
	item.watcher.setFuture(QtConcurrent::run([id, path, token, fetch] {
		Result result;
		try {
			result.original = readCatalog(path, result.existed);
			const auto current = result.original.toStdString();
			const auto update = fetch(id, current, token);
			result.contents = QByteArray::fromStdString(update.contents);
			result.updated = update.updated;
			result.unmapped = update.unmapped;
			result.changed = update.changed;
		} catch (const std::exception& error) {
			result.error = QString::fromUtf8(error.what());
		}
		return result;
	}));
}

void ObservatoryUpdate::check()
{
	checkSource(m_optical);
	checkSource(m_radar);
	emit changed();
}

void ObservatoryUpdate::install(const QString& id)
{
	auto* item = source(id);
	if (!item || !item->available || item->checking) return;
	const auto path = catalogPath(id);
	try {
		bool exists = false;
		const auto current = readCatalog(path, exists);
		if (exists != item->result.existed || current != item->result.original) {
			checkSource(*item);
			item->message = tr("Local data changed; checking again.");
			emit changed();
			return;
		}
		if (!QDir().mkpath(QFileInfo(path).absolutePath()))
			throw std::runtime_error("Cannot create the observatory data directory");
		QSaveFile file(path);
		file.setDirectWriteFallback(false);
		if (!file.open(QIODevice::WriteOnly) || file.write(item->result.contents) != item->result.contents.size() || !file.commit())
			throw std::runtime_error(file.errorString().toStdString());
		item->available = false;
		item->message = tr("Update installed. Restart the application to use it.");
		if (item->result.unmapped)
			item->message += tr(" %1 JPL stations have no verified MPC code.").arg(item->result.unmapped);
		emit changed();
		emit installed();
	} catch (const std::exception& error) {
		item->message = tr("Update failed: %1").arg(QString::fromUtf8(error.what()));
		emit changed();
	}
}
