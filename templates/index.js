function getColorFromValue(value, maxValue = 1552) {
  if (value === -2) return '#000000'; // Negro para el valor especial
  if (value < 0) return '#808080';    // Gris para otros negativos (opcional)

  const hue = (value / maxValue) * 360;
  return `hsl(${hue}, 100%, 50%)`;
}

let puntosActuales = [];
let ultimoMarkerSeleccionado = null;

function buscarCriminalidad() {
    if (!bounds) return alert("Dibuja un área primero");
  
    const data = {
      minLat: bounds.getSouth(),
      minLon: bounds.getWest(),
      maxLat: bounds.getNorth(),
      maxLon: bounds.getEast()
    };
  
    fetch("http://localhost:8080/consultar", {
      method: "POST",
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(data)
    })
    .then(res => res.json())
    .then(puntos => {
        puntosActuales = puntos;
        puntosLayer.clearLayers();
        document.getElementById("rankingPanel").style.display = "none";
  
        puntos.forEach((p) => {
            let color = getColorFromValue(p.color);
            const marker = L.circleMarker([p.lat, p.lon], {
                radius: 4,
                color: color,
                weight: 1,
                fillOpacity: 1,
                fillColor: color,
                stroke: true
            }).addTo(puntosLayer);
    
            p._marker = marker;
    
            marker.on("click", () => {
                if (ultimoMarkerSeleccionado) {
                    ultimoMarkerSeleccionado.setStyle({
                        weight: 1,
                        color: getColorFromValue(ultimoMarkerSeleccionado.options.clusterColor),
                    });
                }
    
                marker.options.clusterColor = color;
    
                marker.setStyle({
                    weight: 3,
                    color: 'black'
                });
    
                ultimoMarkerSeleccionado = marker;
                
                mostrarRankingCompleto(p.indice_h);
            });
        });
    })
    .catch(err => alert("Error en la consulta: " + err));
}
  
function mostrarRankingCompleto(indice_h) {
    const rankingDiv = document.getElementById("rankingLista");
    const panel = document.getElementById("rankingPanel");
    rankingDiv.innerHTML = "";
    panel.style.display = "block";
  
    const ordenados = [...puntosActuales].sort((a, b) => 
      Math.abs(a.indice_h - indice_h) - Math.abs(b.indice_h - indice_h)
    );
  
    ordenados.forEach((p, i) => {
      const div = document.createElement("div");
      div.className = "ranking-item";
  
      const color = getColorFromValue(p.color);
      const bolita = document.createElement("div");
      bolita.className = "bolita";
      bolita.style.backgroundColor = color;
  
      const texto = document.createElement("div");
      texto.innerText = `#${i + 1}: Lat=${p.lat.toFixed(5)}, Lon=${p.lon.toFixed(5)}\nKey=${p.indice_h}`;
  
      div.appendChild(bolita);
      div.appendChild(texto);
      rankingDiv.appendChild(div);
    });
}