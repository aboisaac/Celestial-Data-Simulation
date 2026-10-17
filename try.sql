CREATE DATABASE IF NOT EXISTS astpyth;
USE astpyth;

CREATE TABLE planets (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(50) NOT NULL,
    mass FLOAT,
    distance_from_sun FLOAT,
    moons INT
);

INSERT INTO planets (name, mass, distance_from_sun, moons) VALUES
('Mercury', 0.33, 57.9, 0),
('Venus', 4.87, 108.2, 0),
('Earth', 5.97, 149.6, 1),
('Mars', 0.642, 227.9, 2),
('Jupiter', 1898, 778.6, 79),
('Saturn', 568, 1433.5, 82),
('Uranus', 86.8, 2872.5, 27),
('Neptune', 102, 4495.1, 14);