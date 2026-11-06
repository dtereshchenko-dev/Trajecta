#include "HttpClient.h"
#include <gtest/gtest.h>

TEST(Http, RejectsSchemesAndEmbeddedNull)
{
	for (const auto* url : {"file:///etc/passwd", "ftp://example.org/a", "example.org/a"})
		EXPECT_THROW((void)trajecta::http::get(url, 1, 100, {}), std::invalid_argument);
	EXPECT_THROW((void)trajecta::http::get(std::string("https://example.org\0hidden", 26), 1, 100, {}),
	             std::invalid_argument);
}

TEST(Http, RejectsLimitsBeforeNetwork)
{
	EXPECT_THROW((void)trajecta::http::get("https://example.org", 0, 100, {}), std::invalid_argument);
	EXPECT_THROW((void)trajecta::http::get("https://example.org", 3601, 100, {}), std::invalid_argument);
	EXPECT_THROW((void)trajecta::http::get("https://example.org", 1, 0, {}), std::invalid_argument);
}

TEST(Http, CancelledBeforeNetwork)
{
	std::stop_source stop;
	stop.request_stop();
	EXPECT_THROW((void)trajecta::http::get("https://example.org", 1, 100, stop.get_token()), trajecta::http::Cancelled);
}

TEST(Http, RejectsInvalidRangesBeforeNetwork)
{
	using trajecta::http::ByteRange;
	EXPECT_THROW((void)trajecta::http::get("https://example.org", 1, 100, {}, ByteRange{100, 50}),
	             std::invalid_argument);
	EXPECT_THROW((void)trajecta::http::get("https://example.org", 1, 100, {}, ByteRange{0, 100}),
	             std::invalid_argument);
}
