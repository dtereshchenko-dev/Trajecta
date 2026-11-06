pragma Singleton
import QtQuick

QtObject {
    readonly property var planets: [
        { full: "Mercury", bit: 0,  color: "#a0a0a0" },
        { full: "Venus",   bit: 1,  color: "#e8c46a" },
        { full: "Earth",   bit: 2,  color: "#4a7dff" },
        { full: "Mars",    bit: 3,  color: "#e87461" },
        { full: "Jupiter", bit: 4,  color: "#d4a373" },
        { full: "Saturn",  bit: 5,  color: "#f0d080" },
        { full: "Uranus",  bit: 6,  color: "#7ec8e3" },
        { full: "Neptune", bit: 7,  color: "#5b7dd4" },
        { full: "Moon",    bit: 11, color: "#c8c0b0" },
        { full: "Pluto",   bit: 12, color: "#a9a8a1" }
    ]

    readonly property var asteroids: [
        { full: "Ceres",  bit: 8,  color: "#c8b690" },
        { full: "Pallas", bit: 9,  color: "#a5b4d4" },
        { full: "Vesta",  bit: 10, color: "#d4bb74" }
    ]

    readonly property color sb441Color: "#7ad8c0"
}
