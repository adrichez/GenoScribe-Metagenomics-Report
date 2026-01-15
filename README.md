<p align="center">
  <img src="assets/cover.png" alt="GenoScribe Banner" style="width:100%; border-radius: 10px;">
</p>

<div align="center">
  <h1><span style="color: #6caef1;">Ejemplo de Informe de Metagenómica generado con GenoScribe</span></h1>

  <hr style="border:none; height:0.3px; background-color:#777; width:65%; margin:30px auto 35px auto;">

  <p>
    <a href="https://www.r-project.org/"><img src="https://img.shields.io/badge/R-276DC3?style=flat&logo=r&logoColor=white" alt="R"></a>
    <a href="https://quarto.org/"><img src="https://img.shields.io/badge/Quarto-008080?style=flat&logo=quarto&logoColor=white" alt="Quarto"></a>
    <a href="https://www.python.org/"><img src="https://img.shields.io/badge/Python-3776AB?style=flat&logo=python&logoColor=white" alt="Python"></a>
    <a href="https://www.gnu.org/software/bash/"><img src="https://img.shields.io/badge/Bash-4EAA25?style=flat&logo=gnu-bash&logoColor=white" alt="Bash"></a>
    <a href="https://nextflow.io/"><img src="https://img.shields.io/badge/Nextflow-FF6F00?style=flat&logo=nextflow&logoColor=white" alt="Nextflow"></a>
    <a href="https://www.w3.org/html/"><img src="https://img.shields.io/badge/HTML-E34F26?style=flat&logo=html5&logoColor=white" alt="HTML"></a>
    <a href="https://www.w3.org/Style/CSS/"><img src="https://img.shields.io/badge/CSS-1572B6?style=flat&logo=css3&logoColor=white" alt="CSS"></a>
    <a href="https://www.javascript.com/"><img src="https://img.shields.io/badge/JavaScript-F7DF1E?style=flat&logo=javascript&logoColor=black" alt="JavaScript"></a>
    <a href="https://shiny.posit.co/"><img src="https://img.shields.io/badge/Shiny-FF4400?style=flat&logo=r&logoColor=white" alt="Shiny"></a>
    <a href="https://apptainer.org/"><img src="https://img.shields.io/badge/Apptainer-0052CC?style=flat&logo=apptainer&logoColor=white" alt="Apptainer"></a>
    <a href="https://www.docker.com/"><img src="https://img.shields.io/badge/Docker-2496ED?style=flat&logo=docker&logoColor=white" alt="Docker"></a>
    <a href="https://github.com/"><img src="https://img.shields.io/badge/GitHub-181717?style=flat&logo=github&logoColor=white" alt="GitHub"></a>
  </p>

  <p>
    <a href="## 📄 Descripción">Descripción</a> •
    <a href="## 📂 Estructura del repositorio">Estructura</a> • 
    <a href="## ⚙️ Requisitos">Requisitos</a> • 
    <a href="## 💻 Instalación">Instalación</a> • 
    <a href="## 🚀 Uso">Uso</a> • 
    <a href="## 🎬 Demostración">Demostración</a> • 
    <a href="## 📝 Notas">Notas</a> • 
    <a href="## 🔗 Informe Web">Informe Web</a> •
    <a href="## 📬 Contacto">Contacto</a>
  </p>
</div>








<br>

---

## 📄 Descripción

Este repositorio es un **ejemplo de informe de Metagenómica generado con GenoScribe**, diseñado para mostrar de manera accesible las capacidades del sistema de generación de informes.  

El objetivo principal es **servir como demostración**, alojando el informe en una página web para que pueda consultarse desde un enlace, mostrando de forma interactiva los resultados de un análisis metagenómico.  

La **Metagenómica** es una técnica de secuenciación de **ADN microbiano** que permite estudiar la **composición, diversidad y potencial funcional** de comunidades microbianas en una muestra biológica o ambiental, sin necesidad de aislar y cultivar cada organismo individualmente.  

🔹 **Aplicaciones principales de la Metagenómica:**  
- 📊 Comparación de comunidades microbianas entre condiciones diferentes (ej. microbiota de individuos sanos vs. enfermos, suelos tratados vs. no tratados).  
- 🔍 Detección de cambios en abundancia de microorganismos específicos (bacterias, arqueas, virus, hongos).  
- 🌱 Descubrimiento de rutas metabólicas y funciones microbianas activas o potenciales en un ecosistema o en el cuerpo humano.  
- 🧪 Evaluación de diversidad microbiana y ecología microbiana, comprendiendo interacciones entre microorganismos y su entorno.  
- 💡 Identificación de genes de interés biotecnológico o clínico, como genes de resistencia a antibióticos o enzimas industriales.  

⚡️ **Ventajas de GenoScribe:**  
- Genera **informes interactivos** con gráficos, tablas y resúmenes claros.  
- Permite **reproducibilidad**: cualquier investigador puede ejecutar el análisis con los mismos parámetros y obtener el mismo informe.  
- Facilita que la **bioinformática sea más accesible**, transformando datos complejos en **información visual y compartible**.  

> 💡 **En resumen:** GenoScribe no solo ejecuta análisis, sino que **traduce la complejidad en información útil y comunicable**.








<br>

---

## 📂 Estructura del repositorio

```plaintext
GenoScribe-Metagenomics-Report
├── assets
│   ├── cover.png
│   └── demostration.gif
├── docs
│   ├── resources
│   ├── site_libs
│   └── index.html
└── README.md
```

**Explicación:**

* `assets/` contiene imágenes y GIFs para la portada y demostración.
* `docs/` aloja el informe web generado (`index.html`) y sus recursos (`resources/` y `site_libs/`).
* `README.md` documenta el proyecto.

El informe se genera como un **sitio web auto-contenido**, que puede visualizarse y compartirse de manera independiente:

```plaintext
report/
├── index.html  # Página principal del informe
├── resources/  # Imágenes, scripts, estilos
└── site_libs/  # Librerías necesarias para interactividad
```

Cada informe incorpora un **Mini Chatbot RAG** para interactuar con la información contenida y obtener datos adicionales sobre GenoScribe y análisis bioinformáticos.
Si se quisiera ejecutar localmente, se usaría un script en `resources/1-essential/3-scripts/5-python/run_report_server.sh` para levantar un servidor local y abrir el informe con todas las funcionalidades.








<br>

---

## 💡 Demostración

El pipeline de **Metagenómica** en GenoScribe genera un **informe interactivo, reproducible y auto-contenido**, estructurado en pestañas que permiten explorar los resultados de manera progresiva y clara.

<p align="center">
  <img src="assets/demostration.gif" alt="GenoScribe Demonstration" width="65%" style="border-radius: 10px;">
</p>

<br>

Este GIF muestra la **estructura y funcionalidad del informe**, incluyendo gráficos, tablas y resúmenes de composición y funcionalidad microbiana.






<br>

---

## 🎬 Proceso de Generación del Informe Interactivo

Si desea ver **todo el proceso real** de creación de un informe interactivo de Bulk RNA-Seq con GenoScribe —desde la descarga del repositorio hasta la generación final del reporte— puede consultar el siguiente video demostrativo:

<p align="center">
  <a href="https://youtu.be/ddT5yJihnVE" target="_blank">
    <img src="https://img.youtube.com/vi/ddT5yJihnVE/maxresdefault.jpg" alt="Video Demostrativo de Generación del Informe" width="65%" style="border-radius: 10px;">
  </a>
</p>

<p align="center">
  👉 <a href="https://youtu.be/ddT5yJihnVE" target="_blank"><strong>Ver video en YouTube</strong></a>
</p>






<br>

---

## 🔗 Informe Web

Consulta la versión web interactiva del informe de ejemplo: 👉 [https://adrichez.github.io/GenoScribe-Metagenomics-Report/](https://adrichez.github.io/GenoScribe-Metagenomics-Report/)  

**Objetivo principal:** que el informe esté disponible en la web para su **consulta y demostración**.






<br>

---

## 📝 Notas

* Este repositorio es **principalmente demostrativo**, para mostrar cómo se ve un informe de GenoScribe para Metagenómica.  
* Para más información sobre cómo se generó este informe, los pasos, parámetros y herramientas involucradas, se puede consultar la **Guía de Usuario de GenoScribe para Metagenómica**, disponible en el repositorio principal: [https://github.com/adrichez/GenoScribe](https://github.com/adrichez/GenoScribe).  




<br>

---

## 📬 Contacto

Para consultas, sugerencias o reportes de errores relacionados con este proyecto, puedes contactarnos a través de:

<ul style="list-style:none; padding:0; margin:0; gap:10px;">
  <li style="display:flex; align-items:center; gap:10px; margin-bottom:8px;">
    <a href="mailto:asanca33@gmail.com" style="text-decoration:none;">
      <span style="background-color:#ff6f61;color:white;padding:4px 10px;border-radius:12px;font-weight:bold;">📧 Email</span>
    </a>
    <span><strong>&rArr; Correo electrónico para consultas:</strong> <a href="mailto:asanca33@gmail.com">asanca33@gmail.com</a></span>
  </li>

  <li style="display:flex; align-items:center; gap:10px; margin-bottom:8px;">
    <a href="https://github.com/adrichez/GenoScribe" style="text-decoration:none;">
      <span style="background-color:#6cc644;color:white;padding:4px 10px;border-radius:12px;font-weight:bold;">🐙 GitHub</span>
    </a>
    <span><strong>&rArr; Repositorio y contribuciones:</strong> <a href="https://github.com/adrichez/GenoScribe">https://github.com/adrichez/GenoScribe</a></span>
  </li>

  <li style="display:flex; align-items:center; gap:10px; margin-bottom:8px;">
    <a href="6-info/" style="text-decoration:none;">
      <span style="background-color:#1da1f2;color:white;padding:4px 10px;border-radius:12px;font-weight:bold;">📚 Docs</span>
    </a>
    <span><strong>&rArr; Documentación adicional:</strong> <a href="../6-info/">6-info/</a></span>
  </li>
</ul>

<br>

<hr style="border:none; height:0.3px; background-color:#777; width:65%; margin:20px auto 25px auto;">

<p align="center" style="margin-top:0.5rem; color:#777; font-size:0.95rem;">
  💡 Tus comentarios y contribuciones ayudan a que GenoScribe sea más <strong>robusto, reproducible y accesible</strong> para la comunidad bioinformática.
</p>
