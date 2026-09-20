# Projekt: Numerische Simulation der Kármánschen Wirbelstraße
 
## 1. Projektbeschreibung
Dieses Projekt befasst sich mit der numerischen Simulation und Visualisierung der **Kármánschen Wirbelstraße**. Dabei handelt es sich um ein charakteristisches Phänomen der Strömungsmechanik, bei dem sich hinter einem umströmten Körper periodisch Wirbelablösungen bilden.
 
Die Implementierung erfolgt vollständig in **Python**. Um die physikalischen Gegebenheiten mit hoher numerischer Auflösung abzubilden, wird das Berechnungsgitter in **logarithmischen Polarkoordinaten** formuliert.
 
## 2. Hyperparameter und Simulationsparameter
 
* **Reynolds Zahl ($Re$):** `100` bis `150` (Charakterisiert das Strömungsregime)
* **Gitterauflösung ($N_\xi \times N_\theta$):** `160x160`, `256x256` (Anzahl der Gitterpunkte in radialer und tangentialer Richtung)
* **Zeitschrittweite ($\Delta t$):** `Δt = 0.033` für `m = n = 160`, `Δt = 0.012` für `m = n = 256` (Diskretisierung der Zeitachse)
* **Simulationsdauer ($T_{\text{max}}$):** `180` (dimensionslose Zeiteinheiten)
 Die Gitterweite ergibt sich aus der θ-Diskretisierung zu $h = 2\pi/(m-2)$ und wird in beiden Richtungen verwendet, sodass die Gitterzellen quadratisch sind. Daraus folgt $\xi_{\text{max}} = (n-1)\,h$; bei `n = m = 160` entspricht das einem Rechengebietsradius von $e^{\xi_{\text{max}}} \approx 558$ Zylinderradien.
 
Die angegebenen Zeitschrittweiten liegen nahe an der Stabilitätsgrenze des expliziten Verfahrens. Größere Werte (z. B. `Δt = 0.038` bei `n = m = 160`) führen zur Divergenz der Rechnung.
 
## 3. Algorithmen und Numerische Verfahren
 
### 3.1 Lösung der Poisson-Gleichung (Stromfunktion)
Zur Berechnung der Stromfunktion $\psi$ aus der Wirbelstärke $\omega$ wird in jedem Zeitschritt die Poisson-Gleichung $\nabla^2\psi = -\omega$ gelöst:
 
- **Verfahren:** Direktes Lösen des linearen Gleichungssystems $A\psi = b$. Die Systemmatrix $A$ wird als dünnbesetzte Matrix aufgestellt und einmalig LU-zerlegt, sodass pro Zeitschritt nur noch Vorwärts- und Rückwärtseinsetzen nötig ist.
- **Diskretisierung:** Zentrale Differenzen 2. Ordnung, natürliche Nummerierung $k = i + (j-1)n$.
- **Randbedingungen:**
  - Zylinderrand ($\xi = 0$): $\psi_{1,j} = 0$
  - Fernfeld ($\xi = \xi_{\text{max}}$): $\psi_{n,j} = e^{\xi_n}\sin\theta_j$
  - θ-Richtung: periodisch über zwei Geisterpunkte, $\psi_{i,1} = \psi_{i,m-1}$ und $\psi_{i,m} = \psi_{i,2}$
 ### 3.2 Zeitschrittverfahren der Wirbeltransportgleichung
Die zeitliche Entwicklung des Wirbelfeldes $\omega$ wird über die Wirbeltransportgleichung modelliert:
 
- **Zeitdiskretisierung:** Runge-Kutta-Verfahren 4. Ordnung für die Wirbeltransportgleichung. Die Stromfunktion wird einmal pro Zeitschritt aus der Poisson-Gleichung bestimmt und über alle vier RK-Stufen konstant gehalten (verzögerte ψ-ω-Kopplung). Das Gesamtverfahren ist dadurch erster Ordnung in der Zeit.
- **Diffusionsterm:** Zentrale Differenzen 2. Ordnung
- **Konvektionsterm:** Zentrale Differenzen 2. Ordnung
- **Randbedingungen der Wirbelstärke:**
  - Zylinderrand: $\omega_{1,j} = (\psi_{3,j} - 8\psi_{2,j}) / 2h^2$
  - Fernfeld: $\omega_{n,j} = 0$
  - θ-Richtung: periodisch, analog zur Stromfunktion
 ### 3.3 Anfangsbedingung
Die Simulation startet mit $\omega = 0$ zuzüglich einer kleinen Störung
 
$$\omega_0 = 0.05 \cdot (\cos\theta + \sin\theta) \cdot e^{-2\xi}.$$
 
Diese bricht die Symmetrie zur x-Achse und beschleunigt das Einsetzen der Instabilität erheblich. Bei `Re = 100` und `n = m = 160` wächst die Schwingung ab etwa $t \approx 40$ an und ist ab $t \approx 100$ voll ausgebildet. Ohne Störung entsteht die Wirbelstraße erst bei $t \approx 1100$, was für die hier gewählte Simulationsdauer nicht praktikabel wäre.
 
## 4. Validierung
Zur Validierung wird die Strouhal-Zahl $St = f \cdot D / U$ berechnet. Dazu wird die Wirbelstärke an einem festen Punkt im Nachlauf ($r = 3$, $\theta = 0$) über die Zeit aufgezeichnet und die Ablösefrequenz $f$ aus dem Maximum des Amplitudenspektrums (FFT) der zweiten Hälfte des Zeitsignals bestimmt. Wegen der Entdimensionalisierung mit Zylinderradius und Anströmgeschwindigkeit gilt $St = 2f$.
 
Als Referenz dient die empirische Korrelation nach Roshko:
 
$$St = 0.198 \left(1 - \frac{19.7}{Re}\right)$$
 
Die Auswertung erfolgt im Notebook `notebooks/strouhal_number.ipynb`.
 
## Projektstruktur
```text
Numerik Projekt Karmansche Wirbelstrasse/
├── notebooks/
│   ├── strouhal_number.ipynb
│   └── script_for_visualizing_vorticity_with_scale.py
├── src/
│   ├── grid.py
│   ├── simulation.py
│   └── solver.py
├── readme.md
└── main.py
```
 
## Benutzung
Das Programm kann über die `main.py` ausgeführt werden. Dort lassen sich die Parameter anpassen und die Simulation starten.
 
Im Notebook `notebooks/strouhal_number.ipynb` kann die Simulation ebenfalls gestartet und die Strouhal-Zahl zur Validierung der Simulation berechnet werden.
 
## Abhängigkeiten
- **numpy**: vektorielles Rechnen und sonstige Mathematik
- **matplotlib**: Visualisierung
- **scipy**: Sparse Matrizen und Lösen von LGS
 ## Grundlage
Die numerische Formulierung folgt den Vorlesungsnotizen *Flow Around a Cylinder* von Jeffrey R. Chasnov (HKUST), insbesondere den Kapiteln zu logarithmischen Polarkoordinaten (Lecture 8) sowie zur instationären Strömung (Lectures 16–21).
 
## Autoren
Louis B., Linus D., Gustav H., Maximilian M.
