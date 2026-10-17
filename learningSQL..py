import sqlite3

# 1. Koble til (lager filen planets.db hvis den ikke finnes)
conn = sqlite3.connect("planets.db")
cursor = conn.cursor()

# 2. Lag en tabell
cursor.execute("""
CREATE TABLE IF NOT EXISTS planets (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    mass REAL,
    distance_from_sun REAL,
    moons INTEGER
)
""")

# 3. Sett inn data
planets_data = [
    ("Mercury", 0.33, 57.9, 0),
    ("Venus", 4.87, 108.2, 0),
    ("Earth", 5.97, 149.6, 1),
    ("Mars", 0.642, 227.9, 2),
    ("Jupiter", 1898, 778.6, 79),
    ("Saturn", 568, 1433.5, 82),
    ("Uranus", 86.8, 2872.5, 27),
    ("Neptune", 102, 4495.1, 14)
]

cursor.executemany("INSERT INTO planets (name, mass, distance_from_sun, moons) VALUES (?, ?, ?, ?)", planets_data)

# 4. Lagre endringene
conn.commit()

# 5. Hent og vis data
print("🌍 Alle planeter:")
print("-" * 50)

cursor.execute("SELECT name, mass, moons FROM planets ORDER BY distance_from_sun")
for row in cursor.fetchall():
    name, mass, moons = row
    print(f"{name:10} | Masse: {mass:8} | Måner: {moons}")

# 6. Spørring: Hvilke planeter har mer enn 10 måner?
print("\n🌕 Planeter med mer enn 10 måner:")
cursor.execute("SELECT name, moons FROM planets WHERE moons > 10")
for row in cursor.fetchall():
    print(f"{row[0]} har {row[1]} måner")

# 7. Lukk tilkoblingen
conn.close()