#pragma once

#ifdef _WIN32
#ifdef ORBITAL_CORE_EXPORTS
#define ORBITAL_CORE_API __declspec(dllexport)
#else
#define ORBITAL_CORE_API __declspec(dllimport)
#endif
#else
#define ORBITAL_CORE_API
#endif
