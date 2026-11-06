#pragma once

#include <QJsonObject>
#include <QString>
#include <optional>
#include <stop_token>
#include <vector>
#include <strewn_core/StrewnJob.h>

namespace atmosphere {

	struct EntryEvent {
		double startMjdTdb{};
		double endMjdTdb{};
		std::vector<strewn::EntryState> entries;
	};

	struct EntryArchive {
		QString path;
		QByteArray hash;
		QString photometryPath;
		QByteArray photometryHash;
		std::optional<double> magnitude;
		std::vector<EntryEvent> events;
		QJsonObject provenance;
	};

	EntryArchive readArchive(const QString& path, std::stop_token stop = {});
	void verifyArchive(const EntryArchive& archive);
	void validateEntry(const strewn::EntryState& entry);
	void saveEntries(const QString& path, const std::vector<strewn::EntryState>& entries);
	void saveJson(const QString& path, const QJsonObject& value);
	double parseUtc(const QString& text);
	QString formatTdb(double mjdTdb);

}  // namespace atmosphere
