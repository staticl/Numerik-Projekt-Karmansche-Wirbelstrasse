# Projekt: Numerische Simulation der Kármánschen Wirbelstraße

## 1. Projektbeschreibung
Dieses Projekt befasst sich mit der numerischen Simulation und Visualisierung der **Kármánschen Wirbelstraße**. Dieser Wirbel ist ein charakteristisches Phänomens der Strömungsmechanik, bei dem sich hinter einem umströmten Körper periodisch Wirbelablösungen bilden.

Die Implementierung erfolgt vollständig in **Python**. Um die physikalischen Gegebenheiten mit hoher numerischer Auflösung abzubilden, wird das Berechnungsgitter in **logarithmischen Polarkoordinaten** formuliert.

## 2. Hyperparameter und Simulationsparameter

* **Reynolds Zahl ($Re$):** `100` bis `150` (Charakterisiert das Strömungsregime)
* **Gitterauflösung ($N_\xi \times N_\theta$):** `160x160`, `256x256` (Anzahl der Gitterpunkte in radialer und tangentialer Richtung)
* **Zeitschrittweite ($\Delta t$):** `Δt = 0.038` für `m = n = 160`, `Δt = 0.012` für `m = n = 256` (Diskretisierung der Zeitachse)
* **Simulationsdauer ($T_{\text{max}}$):** `180s` (Gesamte Simulationszeit)

## 3. Algorithmen und Numerische Verfahren

### 3.1 Lösung der Poisson-Gleichung (Stromfunktion)
Zur Berechnung der Stromfunktion $\psi$ aus der Wirbelstärke $\omega$ wird im ersten Schritt die Poisson-Gleichung gelöst:

- **Verfahren:** Lösen eines linearen Gleichungssystems
- **Randbedingungen (periodisch):** Für die Simulation werden hierbei die klassisch periodischen Randbedingungen der Kármánschen Wirbelstraße verwendet.

### 3.2 Zeitschrittverfahren der Wirbeltransportgleichung
Die zeitliche Entwicklung des Wirbelfeldes $\omega$ wird über die Wirbeltransportgleichung modelliert:

- **Zeitdiskretisierung:** Runge-Kutta-Verfahren 4. Ordnung
- **Diffusionsterm:** Zentrale Differenzen 2. Ordnung
- **Konvektion:** Upwind-Schema 3. Ordnung

## Projektstruktur
```text
Numerik Projekt Karmansche Wirbelstrasse/
├── notebooks/
│   └── strouhal_number.ipynb
├── src/
│   ├── grid.py
│   ├── simulation.py
│   └── solver.py
├── README.md
└── main.py
```

## Benutzung
Das Programm kann über die `main.py` ausgeführt werden. Dort lassen sich die Parameter anpassen und die Simulation starten.

Im Notebook `notebooks/strouhal_number.ipynb` kann die Simulation ebenfalls gestartet und die Strouhal-Zahl zur Validierung der Simulation berechnet werden.

## Abhängigkeiten
- **numpy**: vektorielles Rechnen und sonstige Mathematik
- **matplotlib**: Visualisierung
- **scipy**: Sparse Matrizen und Lösen von LGS

## Autoren
Louis Bösenberg, Linus Demuth, Gustav-Theodor Henschel, Maximillian Meran
