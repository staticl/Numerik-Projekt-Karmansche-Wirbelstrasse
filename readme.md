# Projekt: Numerische Simulation der Kármánschen Wirbelstraße

## 1. Projektbeschreibung
Dieses Projekt befasst sich mit der numerischen Simulation und Visualisierung der **Kármánschen Wirbelstraße**. Die **Kármánschen Wirbelstraße** ist ein charakteristisches Phänomens der Strömungsmechanik, bei dem sich hinter einem umströmten Körper periodisch Wirbelablösungen bilden.

Die Implementierung erfolgt vollständig in **Python**. Um die physikalischen Gegebenheiten mit hoher numerischer Auflösung abzubilden, wird das Berechnungsgitter in **logarithmischen Polarkoordinaten** formuliert.

## 2. Hyperparameter und Simulationsparameter

* **Reynolds Zahl ($Re$):** `[die verwendete Reynolds Zahl]` (Charakterisiert das Strömungsregime)
* **Gitterauflösung ($N_\xi \times N_\theta$):** `[Radial x Tangential]` (Anzahl der Gitterpunkte in radiale und tangentiale Richtung)
* **Zeitschrittweite ($\Delta t$):** `[DeltaT]` (Diskretisierung der Zeitachse)
* **Simulationsdauer ($T_{\text{max}}$):** `[t]` (Gesamte Simulationszeit)

## 3. Algorithmen und Numerische Verfahren

### 3.1 Lösung der Poisson Gleichung (Stromfunktion)
Zur Berechnung der Stromfunktion $\psi$ aus der Wirbelstärke $\omega$ wird im ersten Schritt die Poisson Gleichung gelöst:

* **Verfahren:** `Successive Over Relaxation (SOR)`
* **Randbedingungen (bisherig):** `[Hier kurze Beschreibung der bisherigen Randbedingungen]`
* **Randbedingungen (periodisch):** `[Hier kurze Beschreibung der periodischen Randbedingungen]`

### 3.2 Zeitschrittverfahren der Wirbeltransportgleichung
Die zeitliche Entwicklung des Wirbelfeldes $\omega$ wird über die Wirbeltransportgleichung modelliert:

* **Zeitdiskretisierung:** `bisher expliziter Euler (1.Ordnung) -> ändern zu Runge-Kutta 4.Ordnung`
* **Diffusionsterm:** `Zentrale Differenzen 2.Ordnung`
* **Advektions:** `Zentrale Differenzen 2.Ordnung -> ändern zu Upwind-Shema (Unterscheidung Vor- und Rückwärtsdifferenz)`