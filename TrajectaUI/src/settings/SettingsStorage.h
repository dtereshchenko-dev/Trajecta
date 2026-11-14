#pragma once

#include <QString>

class QSettings;

namespace SettingsStorage {
	bool migrate(QSettings& source, const QString& destination, bool removeSource = false,
	             bool* createdEmpty = nullptr);
	bool initialize(bool* createdEmpty = nullptr);
}  // namespace SettingsStorage
