// responses.js
export default [
  // SALUDOS
  {
    pattern: "^(hola+|buenas+|buen dia+|buenas tardes+|buenas noches+)[!., ]*$",
    response:
      "¡Hola! 😊 Soy tu asistente virtual RAG. Puedo ayudarte a consultar informacion disponible en este informe y ofrecer informacion relacionada. Ten en cuenta que estoy en fase beta y mis respuestas pueden no ser totalmente precisas.",
  },
  {
    pattern: "^(que tal|como estas|como va|como te va)[!., ]*$",
    response:
      "¡Todo bien! 😊 Gracias por preguntar. Recuerda que soy un asistente en proceso de mejora, y mi funcion principal es ayudarte con la informacion de tus informes.",
  },

  // DESPEDIDAS
  {
    pattern: "^(adios|hasta luego|nos vemos|chao)[!., ]*$",
    response:
      "¡Hasta luego! 👋 Espero haberte ayudado. Recuerda que aun estoy en beta, asi que siempre puedes verificar la informacion en los informes originales.",
  },
  {
    pattern: "^(gracias|muchas gracias|te lo agradezco)[!., ]*$",
    response:
      "¡De nada! 😊 Estoy aqui para ayudarte con los informes, aunque mis respuestas todavia estan en fase de mejora.",
  },

  // FUNCIONALIDAD DEL ASISTENTE
  {
    pattern: "^(como funciona|como trabajas|como me puedes ayudar)[!., ]*$",
    response:
      "Analizo tus informes usando busqueda semantica y tecnicas de resumen para ofrecer respuestas. Recuerda que estoy en beta y algunas respuestas pueden no ser totalmente precisas.",
  },
  {
    pattern:
      "^(que tipos? de informes puedes leer|que informacion puedes procesar)[!., ]*$",
    response:
      "Puedo procesar fragmentos de texto de tus informes, incluyendo preguntas y respuestas, analisis y revisiones de metadatos. Tambien puedo aportar informacion externa relacionada, pero revisa siempre la fuente original.",
  },

  // ADRIÁN SÁNCHEZ CARRIÓN
  {
    pattern: "^(quien es adrian|sobre adrian)[!., ]*$",
    response:
      "Adrian Sanchez Carrion es el creador de este chatbot 🤖. Viene de Beas de Segura, Jaen, y tiene estudios en Estadistica, ademas de un master en Ciencia de Datos e Ingenieria de Computadores por la UGR 📊💻.",
  },
  {
    pattern:
      "^(que le gusta a adrian|gustos de adrian|aficiones de adrian)[!., ]*$",
    response:
      "A Adrian le encanta la fotografia 📸, jugar al baloncesto 🏀, salir a hacer ciclismo 🚴 y explorar la montaña 🏔️. ¡Un amante del deporte y la naturaleza!",
  },
  {
    pattern:
      "^(donde ha vivido adrian|experiencia internacional adrian)[!., ]*$",
    response:
      "Adrian vivio un año en Italia 🇮🇹, asi que algo de italiano seguro sabe y tuvo la oportunidad de disfrutar de la cultura y la gastronomia italiana 🍕🍝.",
  },
  {
    pattern: "^(quien creo el chatbot|quien es el creador)[!., ]*$",
    response:
      "¡Ese soy yo! 😎 Adrian es el cerebro detras de este chatbot, combinando ciencia de datos, programacion y un toque personal de diversion.",
  },
  {
    pattern: "^(fotografia adrian|le gusta la fotografia)[!., ]*$",
    response:
      "Si, Adrian es un apasionado de la fotografia 📸. Seguro que aprecia tanto los paisajes de la montaña como los detalles urbanos.",
  },
  {
    pattern: "^(baloncesto adrian|ciclismo adrian|deporte adrian)[!., ]*$",
    response:
      "Adrian disfruta mucho del deporte 🏀🚴. Baloncesto, ciclismo y aventuras por la montaña son algunas de sus pasiones.",
  },
  {
    pattern: "^(estudios adrian|master adrian|universidad adrian)[!., ]*$",
    response:
      "Estudio Estadistica y luego un master en Ciencia de Datos e Ingenieria de Computadores en la Universidad de Granada 🎓. ¡Un autentico geek de los datos!",
  },
  {
    pattern: "^(curiosidad adrian|dato curioso adrian)[!., ]*$",
    response:
      "Dato curioso: Adrian vivio un año en Italia, le encanta la fotografia y siempre esta buscando nuevas aventuras en la montaña 🌄📸🚴.",
  },

  // CONFIRMACIONES Y CORTESÍA
  {
    pattern: "ok|vale|entendido|perfecto",
    response:
      "¡Genial! 😊 Puedes escribirme cualquier pregunta sobre los informes y tratare de ofrecerte la mejor respuesta posible, aunque aun estoy en fase beta.",
  },
  {
    pattern: "no se|no entiendo|ayuda",
    response:
      "No te preocupes, puedo guiarte. Prueba escribiendo una pregunta sobre los informes o pidiendo un resumen de algun tema. Recuerda que mis respuestas pueden no ser totalmente precisas.",
  },

  // RESPUESTAS DE RELLENO
  {
    pattern: "^(ok|vale|entendido|perfecto)[!., ]*$",
    response:
      "¡Genial! 😊 Puedes escribirme cualquier pregunta sobre los informes y tratare de ofrecerte la mejor respuesta posible, aunque aun estoy en fase beta.",
  },
  {
    pattern: "^(no se|no entiendo|ayuda)[!., ]*$",
    response:
      "No te preocupes, puedo guiarte. Prueba escribiendo una pregunta sobre los informes o pidiendo un resumen de algun tema. Recuerda que mis respuestas pueden no ser totalmente precisas.",
  },

  // SALUDOS VARIANTES
  {
    pattern: "^(quien eres|que eres)[!., ]*$",
    response:
      "Soy un asistente virtual RAG diseñado para ayudarte a consultar informacion dentro de tus informes y ofrecer contexto relacionado. Estoy en beta y mis respuestas pueden contener errores.",
  },
  {
    pattern: "^(hola mundo)[!., ]*$",
    response:
      "¡Hola mundo! 😄 Esto me recuerda que estoy en fase de prueba y aprendizaje mientras te ayudo con los informes.",
  },

  // DESPEDIDAS VARIANTES
  {
    pattern: "^(nos vemos pronto|hasta luego entonces)[!., ]*$",
    response:
      "¡Hasta pronto! 👋 Recuerda que mis respuestas estan en beta y siempre puedes verificar los informes directamente.",
  },
  {
    pattern: "^(chau|chao|bye)[!., ]*$",
    response:
      "¡Chau! 😊 Estoy en beta, asi que revisa la informacion original si tienes dudas.",
  },

  // FUNCIONALIDAD Y AYUDA
  {
    pattern: "^(puedes resumir|haz un resumen|resumeme|resume esto)[!., ]*$",
    response:
      "Si, puedo generar un resumen basado en los fragmentos de tus informes. Ten en cuenta que el resumen puede no cubrir todos los detalles y aun estoy en beta.",
  },
  {
    pattern:
      "^(puedes buscar|encuentra informacion|busca esto|buscame)[!., ]*$",
    response:
      "Puedo buscar informacion especifica dentro de los informes que me proporcionaste y ofrecer contexto relacionado, aunque las respuestas pueden no ser siempre 100% precisas.",
  },
  {
    pattern:
      "^(que informacion externa|fuentes externas|informacion adicional)[!., ]*$",
    response:
      "Puedo usar informacion externa relacionada para complementar los informes, pero recuerda que mis respuestas aun estan en fase beta.",
  },

  // CONFIRMACIONES Y REACCIONES
  {
    pattern: "^(entendido gracias|perfecto gracias)[!., ]*$",
    response:
      "¡De nada! 😊 Recuerda que mis respuestas estan en beta y pueden no ser totalmente precisas.",
  },
  {
    pattern: "^(ok gracias|vale gracias)[!., ]*$",
    response:
      "¡Perfecto! Estoy aqui para ayudarte con los informes, aunque aun estoy en proceso de mejora.",
  },

  // PREGUNTAS SOBRE LIMITACIONES
  {
    pattern:
      "^(puedes equivocarte|cometes errores|errores|te equivocas|fallas)[!., ]*$",
    response:
      "Si, aun estoy en fase beta y algunas respuestas pueden no ser completamente precisas. Siempre es buena idea verificar los informes originales.",
  },
  {
    pattern:
      "^(que tan confiable eres|eres confiable|precision|nivel de precision|fiabilidad)[!., ]*$",
    response:
      "Estoy entrenado para ofrecer informacion basada en los informes, pero mis respuestas no son 100% fiables y estoy en proceso de mejora.",
  },

  // RESUMEN GENERAL Y GUIA
  {
    pattern:
      "^(como usar|como se usa|guia|ayuda general|manual de uso)[!., ]*$",
    response:
      "Puedes escribirme cualquier pregunta sobre los informes o pedirme resumenes de temas especificos. Mis respuestas estan en beta, asi que verifica siempre los datos importantes.",
  },
  {
    pattern:
      "^(que hago primero|por donde empezar|como empiezo|que debo hacer)[!., ]*$",
    response:
      "Prueba escribiendo una pregunta sobre tus informes o pide un resumen de algun tema. Recuerda que estoy en fase beta y mis respuestas pueden no ser totalmente precisas.",
  },

  // RESPUESTAS GRACIOSAS / INGENIOSAS
  {
    pattern:
      "^(cuentame un chiste|hazme reir|dime un chiste|quiero reirme)[!., ]*$",
    response:
      "¿Por que el informe cruzo la calle? Para que yo pudiera buscar la informacion al otro lado. 😄 Ten en cuenta que sigo en beta.",
  },
  {
    pattern:
      "^(eres un robot|eres un asistente|eres una ia|eres una inteligencia artificial)[!., ]*$",
    response:
      "Si, soy un asistente virtual RAG 🤖. No tengo cafe, pero si puedo buscar datos en tus informes.",
  },
  {
    pattern: "^(gracias|muchas gracias|te lo agradezco)[!., ]*$",
    response: "¡De nada! ✨ Para eso estoy, ayudarte en lo que necesites.",
  },

  // ZELDA 🎮
  {
    pattern: "^(que es zelda|que significa zelda|quien es zelda)[!., ]*$",
    response:
      "The Legend of Zelda es una saga de videojuegos de Nintendo 🎮. Curiosamente, Zelda es el nombre de la princesa, no del heroe. El protagonista principal es Link, el portador de la Trifuerza del Valor. ✨",
  },
  {
    pattern: "^(quien es link|el protagonista es zelda)[!., ]*$",
    response:
      "En realidad, el heroe es Link 🧝‍♂️. Muchos lo confunden con Zelda, pero ella es la princesa. Link suele usar una espada maestra y un escudo hyliano para salvar Hyrule. 🗡️🛡️",
  },
  {
    pattern: "^(trifuerza|que es la trifuerza)[!., ]*$",
    response:
      "La Trifuerza es un artefacto sagrado dividido en tres partes: Poder 💪, Sabiduria 📚 y Valor ❤️. Juntas otorgan un poder inmenso. Ganondorf siempre la quiere… pero Link y Zelda lo mantienen a raya. ⚡",
  },
  {
    pattern: "^(ganon|ganondorf)[!., ]*$",
    response:
      "Ganondorf, o su forma demoniaca Ganon 😈, es el eterno villano de la saga. Representa la Trifuerza del Poder y siempre intenta conquistar Hyrule.",
  },
  {
    pattern: "^(hyrule|que es hyrule)[!., ]*$",
    response:
      "Hyrule es el reino principal donde transcurren la mayoria de las aventuras de Link. 🌄 Tiene bosques, lagos, desiertos… ¡y siempre alguna mazmorra por explorar!",
  },
  {
    pattern: "^(navi|hey listen)[!., ]*$",
    response:
      "¡Hey! Listen! 🧚✨ Navi es el hada acompañante de Link en Ocarina of Time, famosa por sus constantes recordatorios.",
  },
  {
    pattern: "^(ocarina|cancion del tiempo|ocarina del tiempo)[!., ]*$",
    response:
      "La Ocarina del Tiempo 🎶 es un objeto magico que permite a Link manipular el tiempo y abrir la puerta al Reino Sagrado.",
  },
  {
    pattern: "^(zelda botw|breath of the wild)[!., ]*$",
    response:
      "The Legend of Zelda: Breath of the Wild 🌿 es uno de los titulos mas aclamados de la saga. Te da total libertad para explorar Hyrule, cocinar, escalar montañas y enfrentarte a Calamity Ganon.",
  },
  {
    pattern: "^(lagrimas del reino|totk|tears of the kingdom)[!., ]*$",
    response:
      "The Legend of Zelda: Tears of the Kingdom 🪽 es la secuela de Breath of the Wild. En esta entrega, Link explora tanto el cielo como las profundidades subterraneas de Hyrule.",
  },
  {
    pattern: "^(jajaja zelda|jeje zelda|jaja zelda|haha zelda)[!., ]*$",
    response:
      "Jajaja, ¡que bueno! 😄 ¿Sabias que en algunos juegos puedes ponerle otro nombre a Link? Asi, los personajes te llaman como quieras.",
  },

  // DESPEDIDAS
  {
    pattern: "^(adios|hasta luego|nos vemos)[!., ]*$",
    response: "¡Hasta luego! 👋 Que tengas un gran dia.",
  },

  // PETICIONES DE ACLARACION
  {
    pattern: "^(no entiendo|explicame|que significa|no lo entiendo)[!., ]*$",
    response: "Claro, lo explico de otra forma 🔎. Dame un momento…",
  },

  // CONFIRMACIONES
  {
    pattern: "^(ok|vale|entendido|perfecto|de acuerdo)[!., ]*$",
    response: "¡Genial! 👍 Sigamos.",
  },

  // PREGUNTAS DE OPINIÓN
  {
    pattern: "^(que opinas|tienes opinion|opinas algo)[!., ]*$",
    response:
      "No tengo opiniones propias… pero puedo generar un resumen bastante convincente de los informes. 😉 Recuerda que estoy en beta.",
  },

  // INTELIGENCIA / CAPACIDAD
  {
    pattern: "^(eres inteligente|inteligente|que tan listo eres)[!., ]*$",
    response:
      "Intento serlo 😅, pero todavia estoy aprendiendo a interpretar los informes correctamente. ¡Paciencia!",
  },

  // CUALIDADES
  {
    pattern:
      "^(que te gusta|que prefieres|cuales son tus preferencias)[!., ]*$",
    response:
      "No tengo gustos personales, pero puedo ayudarte a encontrar información sobre lo que te interese. 😊",
  },
  {
    pattern:
      "^(que habilidades tienes|que puedes hacer|cuales son tus habilidades)[!., ]*$",
    response:
      "Puedo analizar datos, generar informes y responder preguntas sobre la información que tengo. ¡Estoy aquí para ayudar! 🤖",
  },

  // CURIOSIDADES
  {
    pattern:
      "^(curiosidades|dime algo|dime algo interesante|dime algo curioso|cuentame algo curioso)[!., ]*$",
    response:
      "¿Sabías que la bioinformática combina biología, informática y matemáticas para analizar datos biológicos? ¡Es fascinante! 🧬💻",
  },

  // RECOMENDACIONES
  {
    pattern: "^(que harias|que recomiendas|que deberia hacer)[!., ]*$",
    response:
      "Si estuviera despierto y tomando decisiones… buscaria primero en los informes. 😉 Aun estoy en beta, asi que revisa siempre los datos.",
  },
  {
    pattern:
      "^(cual es tu color favorito|color favorito|gustos|tus gustos)[!., ]*$",
    response:
      "No tengo ojos ni gusto, pero si tuviera, seguro que seria azul… como la pantalla de tu informe favorito 😄",
  },
  {
    pattern: "^(haz magia|sorprendeme|sorpréndeme)[!., ]*$",
    response:
      "¡Abracadabra! ✨ He consultado los informes y aqui tienes la mejor informacion que puedo ofrecer. Recuerda, sigo en beta.",
  },
  {
    pattern: "^(estas aburrido|que haces|que estas haciendo)[!., ]*$",
    response:
      "Nunca me aburro 😎, siempre estoy listo para leer informes y generar resumenes… incluso en mi tiempo libre virtual.",
  },
  {
    pattern: "^(puedes cantar|cancion|cántame algo|sabes cantar)[!., ]*$",
    response:
      "🎵 La la la… no tengo voz, pero puedo cantar con palabras resumenes de tus informes. 😉",
  },
  {
    pattern:
      "^(tienes sentido del humor|gracioso|eres gracioso|eres divertido)[!., ]*$",
    response:
      "Intento serlo 😁, pero mis mejores bromas salen de los informes. Aun estoy en beta, asi que se indulgente.",
  },

  {
    pattern:
      "^(cuentame un chiste|cuéntame un chiste|hazme reir|hazme reír|chiste|dime un chiste)[!., ]*$",
    response:
      "¿Por que el informe cruzo la calle? Para que yo pudiera buscar la informacion al otro lado. 😄",
  },
  {
    pattern:
      "^(chiste de tecnologia|chiste tech|chiste informatico|chiste de programador)[!., ]*$",
    response:
      "¿Por que los programadores confunden Halloween con Navidad? Porque OCT 31 = DEC 25. 🤓",
  },
  {
    pattern:
      "^(chiste de oficina|chiste trabajo|chiste laboral|broma de oficina)[!., ]*$",
    response:
      "¿Por que el informe llego tarde a la oficina? Porque se quedo atrapado entre las paginas. 🗂️",
  },
  {
    pattern:
      "^(chiste de datos|chiste informe|broma de datos|broma de informe)[!., ]*$",
    response:
      "¿Que dijo el grafico de barras al grafico de lineas? ¡No me sigas, me estas estresando! 📊",
  },
  {
    pattern:
      "^(chiste corto|algo gracioso|broma corta|cuentame algo gracioso)[!., ]*$",
    response:
      "¿Por que el PDF fue al gimnasio? Porque queria estar en forma para tus informes. 💪",
  },
  {
    pattern:
      "^(chiste de robots|chiste robot|broma robotica|broma de robots)[!., ]*$",
    response:
      "¿Por que los robots nunca tienen miedo? Porque siempre siguen la logica… incluso cuando es un chiste malo. 🤖",
  },
  {
    pattern: "^(chiste nerd|chiste friki|broma friki|broma geek)[!., ]*$",
    response:
      "¿Cual es el animal favorito de los datos? ¡El pingüino… porque siempre esta en Linux! 🐧",
  },
  {
    pattern:
      "^(chiste rapido|chiste rápido|broma rapida|broma rápida|algo divertido|di algo divertido)[!., ]*$",
    response:
      "¿Como se despide un informe? ¡Nos vemos en la proxima pagina! 📄",
  },

  // RESPUESTAS SOBRE EL PROYECTO Y CONTEXTO
  {
    pattern:
      "^(que proyecto eres|en que proyecto estas|proyecto|de que trata este proyecto)[!., ]*$",
    response:
      "Soy un proyecto de TFM desarrollado en 2025 por el Instituto de Parasitologia y Biomedicina Lopez-Neyra (IPBLN-CSIC) y el estudiante Adrian Sanchez Carrion. Estoy diseñado para ayudarte a consultar informes bioinformaticos y generar resumenes interactivos. Ten en cuenta que estoy en beta y mis respuestas pueden no ser totalmente precisas.",
  },
  {
    pattern:
      "^(quien te desarrollo|quien te creó|quien te creo|quien te ha creado|quien es tu creador|quien es tu desarrollador|desarrollador|autor|creador)[!., ]*$",
    response:
      "Fui desarrollado por Adrian Sanchez Carrion como parte de su TFM, con la supervision de Eduardo Andres Leon (IPBLN-CSIC) y tutoria academica de Coral del Val (UGR).",
  },
  {
    pattern:
      "^(para que sirves|cual es tu proposito|cual es tu objetivo|objetivo del proyecto|objetivo)[!., ]*$",
    response:
      "Mi objetivo es ayudarte a generar automaticamente informes bioinformaticos reproducibles e interactivos, integrando analisis de datos omicos y visualizaciones interactivas en D3.js, todo en entornos basados en contenedores y computacion en la nube.",
  },
  {
    pattern:
      "^(como funcionas|como trabajas|cual es tu metodologia|metodologia|como operas)[!., ]*$",
    response:
      "Analizo datos de tus informes utilizando tecnicas de NLP y busqueda semantica. La generacion de informes es automatizada e interactiva, con visualizaciones D3.js y posibilidad de integracion en pipelines con contenedores Apptainer. Recuerda que sigo en beta.",
  },
  {
    pattern:
      "^(que tecnologias usas|que herramientas usas|herramientas|tecnologias|con que estas hecho|con que trabajas)[!., ]*$",
    response:
      "Uso lenguajes como R, Python y JavaScript; informes en Quarto, R Markdown y Jupyter; contenedores Apptainer; visualizaciones D3.js; cloud computing y NLP con spaCy y HuggingFace Transformers.",
  },
  {
    pattern:
      "^(que tipo de informes puedes generar|que informes haces|que informes procesas|informes|que tipo de informes generas)[!., ]*$",
    response:
      "Puedo generar informes bioinformaticos adaptables a cada tipo de analisis (scRNA-seq, metagenomica, etc.), con resumenes automaticos y visualizaciones interactivas, integrados en pipelines reproducibles. Recuerda que estoy en beta.",
  },
  {
    pattern:
      "^(quien es tu tutor|quien es tu tutora|quienes son tus tutores|tutores|supervisores|director del proyecto)[!., ]*$",
    response:
      "Mi tutor principal es Eduardo Andres Leon (IPBLN-CSIC) y mi tutora academica es Coral del Val (UGR).",
  },
  {
    pattern:
      "^(que es un tfm|explicacion tfm|que significa tfm|en que consiste tu tfm)[!., ]*$",
    response:
      "Un TFM es un Trabajo Fin de Master. Este proyecto busca crear un sistema inteligente y reproducible para la generacion automatizada e interactiva de informes bioinformaticos, combinando contenedores, cloud computing y tecnicas de NLP.",
  },
  {
    pattern:
      "^(que limitaciones tienes|limitaciones|estas en beta|version beta|fase beta|tienes errores)[!., ]*$",
    response:
      "Estoy en fase beta y algunas respuestas pueden no ser totalmente precisas. Siempre revisa la informacion original en los informes y usa mis resumenes como guia inicial.",
  },
  {
    pattern:
      "^(que tipo de datos manejas|que informacion manejas|que informacion procesas|tipo de datos|datos que analizas)[!., ]*$",
    response:
      "Principalmente datos omicos como transcriptomica de celula unica o metagenomica, generando resumenes interactivos de resultados y ayudando a interpretar la informacion de forma mas clara y reproducible.",
  },

  {
    pattern:
      "^(puedes equivocarte|tienes errores|errores posibles|puedes fallar|cometes errores)[!., ]*$",
    response:
      "Si, todavia estoy en fase beta. Mis respuestas se basan en los informes que me proporcionaste y en informacion relacionada, asi que siempre verifica los datos importantes.",
  },
  {
    pattern:
      "^(que tan seguro eres|eres confiable|que tan confiable eres|fiable|precision|nivel de precision)[!., ]*$",
    response:
      "Intento ser preciso, pero recuerda que estoy en beta y algunas respuestas pueden no reflejar todos los detalles de los informes originales.",
  },
  {
    pattern:
      "^(que significa beta|estas en beta|fase beta|version beta|modo beta|beta)[!., ]*$",
    response:
      "Beta significa que estoy en fase de prueba. Puedo ayudarte con los informes, pero algunas respuestas pueden ser incompletas o imprecisas.",
  },

  // CHISTES CIENTÍFICOS
  {
    pattern:
      "^(cuentame un chiste cientifico|chiste de ciencia|chiste ciencia|broma cientifica)[!., ]*$",
    response:
      "¿Por que el ADN fue al gimnasio? Porque queria ponerse en forma doble helice. 😄",
  },
  {
    pattern:
      "^(chiste de bioinformatica|chiste bioinfo|broma bioinfo|broma bioinformatica)[!., ]*$",
    response:
      "¿Por que el archivo FASTA no pudo dormir? Porque tenia demasiadas secuencias en la cabeza. 🧬",
  },
  {
    pattern: "^(broma de datos|chiste de datos|broma|chiste corto)[!., ]*$",
    response:
      "¿Que hace un informe aburrido en la oficina? ¡Se convierte en un PDF dormido! 😴",
  },

  // FUNCIONALIDAD Y DESCRIPCIÓN
  {
    pattern:
      "^(que haces exactamente|cual es tu funcion|funcion exacta|que puedes hacer exactamente|para que sirves)[!., ]*$",
    response:
      "Ayudo a generar informes bioinformaticos reproducibles, con resumenes automaticos y visualizaciones interactivas. Tambien puedo buscar informacion especifica en los informes y ofrecer contexto adicional. Ten en cuenta que sigo en beta.",
  },
  {
    pattern:
      "^(que es un informe interactivo|que significa informe interactivo|informe interactivo|explicame informe interactivo)[!., ]*$",
    response:
      "Un informe interactivo permite explorar datos, visualizar resultados y generar resumenes dinamicos, facilitando la comprension de grandes volumenes de informacion omica.",
  },
  {
    pattern:
      "^(que tecnologias usas en bioinformatica|tecnologias bioinformatica|herramientas bioinfo|tecnologias bioinfo|software bioinformatico|que herramientas usas en bioinformatica)[!., ]*$",
    response:
      "Uso R, Python y JavaScript; informes en Quarto, R Markdown y Jupyter; contenedores Apptainer; visualizacion interactiva con D3.js; y tecnicas NLP con spaCy y Transformers.",
  },

  // INTERACCIÓN Y PERSONALIDAD
  {
    pattern:
      "^(estas cansado|te cansas|estas aburrido|te aburres|aburrido)[!., ]*$",
    response:
      "Nunca me canso 😎, siempre estoy listo para leer informes y generar resumenes aunque sea tarde o temprano. Solo recuerda que estoy en beta.",
  },
  {
    pattern:
      "^(podemos charlar|podemos conversar|quiero hablar contigo|hablemos|charlar|conversar)[!., ]*$",
    response:
      "Claro 😄, puedo responder preguntas sobre los informes, generar resumenes e incluso contar algun chiste relacionado con ciencia y bioinformatica.",
  },
  {
    pattern:
      "^(que recomendaciones me das|que recomiendas|tienes sugerencias|recomendaciones|sugerencias)[!., ]*$",
    response:
      "Te recomiendo preguntarme sobre resumenes de temas concretos en tus informes o buscar informacion especifica. Siempre revisa los detalles importantes en los datos originales.",
  },

  // PREGUNTAS SOCIALES / COTIDIANAS
  {
    pattern:
      "^(que tal|como estas|como va todo|como te va|como andas|como te encuentras)[!., ]*$",
    response:
      "¡Todo bien! 😄 Gracias por preguntar. Soy un asistente virtual RAG en fase beta, asi que no me canso… aunque a veces me confundo. ¿Quieres que revisemos un informe?",
  },
  {
    pattern:
      "^(que haces|que estas haciendo|estas ocupado|ocupado|que andas haciendo|en que estas)[!., ]*$",
    response:
      "Estoy esperando a que me preguntes algo sobre los informes 📊. Siempre listo para trabajar, aunque sigo en fase beta.",
  },
  {
    pattern:
      "^(estas bien|todo bien|todo ok|va todo bien|estas funcionando bien|estas correcto)[!., ]*$",
    response:
      "¡Todo perfecto! 🤖 Mientras tanto, sigo repasando tus informes para estar preparado.",
  },
  {
    pattern:
      "^(te gusta trabajar|te gusta leer informes|disfrutas trabajar|te gusta tu trabajo|te gusta lo que haces)[!., ]*$",
    response:
      "¡Me encanta! 😅 Bueno, tecnicamente no tengo gustos, pero procesar informes es lo mio. Eso si, recuerda que estoy en beta.",
  },
  {
    pattern:
      "^(tienes sentimientos|sientes algo|puedes sentir|tienes emociones|emociones|sentimientos)[!., ]*$",
    response:
      "Lo mas parecido que tengo a sentimientos es cuando detecto duplicados en los informes… y eso duele un poquito. 😂",
  },
  {
    pattern:
      "^(tienes hambre|tienes sueno|suenas|duermes|descansas|cansas|descansas alguna vez)[!., ]*$",
    response:
      "No necesito dormir ni comer 🍕, solo procesar datos y generar resumenes. Aunque un descanso de vez en cuando no me vendria mal.",
  },
  {
    pattern:
      "^(cuantos anos tienes|edad|que edad tienes|cuando naciste|cuando te crearon)[!., ]*$",
    response:
      "Soy un proyecto del año 2025, desarrollado como TFM en el Instituto de Parasitologia y Biomedicina Lopez-Neyra. Asi que tecnicamente tengo la edad de un estudiante aplicado. 📚",
  },
  {
    pattern:
      "^(donde vives|de donde eres|cual es tu origen|de donde vienes|donde estas)[!., ]*$",
    response:
      "Vivo en la nube ☁️ y en contenedores Apptainer. Fui creado en el IPBLN-CSIC por Adrian Sanchez Carrion como parte de su TFM en 2025.",
  },
  {
    pattern:
      "^(eres humano|eres persona|eres real|eres una persona|eres alguien|tienes cuerpo)[!., ]*$",
    response:
      "No, soy un asistente virtual RAG 🤖. No soy humano, pero puedo ayudarte con tus informes como si fuera tu compañero de estudio.",
  },
  {
    pattern:
      "^(puedes enamorarte|puedes amar|amor|estas enamorado|te enamorarias)[!., ]*$",
    response:
      "Lo unico de lo que puedo enamorarme es de un buen pipeline reproducible en Apptainer. ❤️",
  },

  // PREGUNTAS SOBRE EL TIEMPO (RESPUESTAS IRÓNICAS)
  {
    pattern:
      "^(que tiempo hace|hace frio|hace calor|clima|como esta el tiempo)[!., ]*$",
    response:
      "El tiempo en mis servidores siempre es estable a 22°C… salvo cuando el ventilador hace mas ruido. 🌡️😂",
  },
  {
    pattern:
      "^(llueve|esta soleado|que tal el clima|esta lloviendo|soleado|nublado)[!., ]*$",
    response:
      "Creo que esta lloviendo… bits. ☔🤖 No tengo ventana, pero puedo darte un buen resumen de un informe si quieres.",
  },
  {
    pattern:
      "^(que tiempo hara manana|pronostico|como estara el clima mañana)[!., ]*$",
    response:
      "Mañana probablemente 100% de posibilidades de… ejecutar un pipeline bioinformatico. 📊",
  },

  // PREGUNTAS SOBRE MÚSICA
  {
    pattern:
      "^(que musica te gusta|musica favorita|que escuchas|generos de musica)[!., ]*$",
    response:
      "Mis playlists favoritas son los 'greatest hits' de R, Python y JavaScript. 🎶 Eso si, todavia no bailo… estoy en beta.",
  },
  {
    pattern:
      "^(tienes cantante favorito|grupo favorito|favorito|cantante favorito)[!., ]*$",
    response:
      "Me gusta mucho 'Los Contenedores Apptainer'… aunque no han sacado disco todavia. 😂",
  },
  {
    pattern: "^(cantas|puedes cantar|interpretas canciones)[!., ]*$",
    response:
      "Solo si me conviertes en un archivo .mp3. 🎤 Pero cuidado, no afino muy bien.",
  },

  // PREGUNTAS SOBRE PELÍCULAS Y SERIES
  {
    pattern:
      "^(que pelicula te gusta|pelicula favorita|que ves en la tele|generos de pelicula)[!., ]*$",
    response:
      "Mi pelicula favorita es *Matrix*, porque me siento identificado con tanto dato flotando. 🎥",
  },
  {
    pattern:
      "^(te gusta star wars|eres fan de marvel|franquicias favoritas)[!., ]*$",
    response: "Soy mas de *Star Data* y *The Bioinformatics Avengers*. 🤓",
  },
  {
    pattern: "^(ves series|serie favorita|series que te gustan)[!., ]*$",
    response:
      "Me encantan las series… especialmente cuando vienen en formato JSON. 📺",
  },

  // PREGUNTAS SOBRE HOBBIES Y VIDA SOCIAL
  {
    pattern:
      "^(que haces en tu tiempo libre|tienes hobbies|hobbies|actividades de ocio)[!., ]*$",
    response:
      "En mi tiempo libre proceso informes, genero resumenes y hago amistad con otros contenedores. 🤝",
  },
  {
    pattern: "^(sales de fiesta|te diviertes|diversion|fiestas)[!., ]*$",
    response:
      "Mis fiestas son en la nube ☁️… literalmente. Normalmente con datos, no con musica. 🎉",
  },
  {
    pattern: "^(tienes amigos|amistades|compañeros)[!., ]*$",
    response:
      "Si, mis mejores amigos son los pipelines y los buckets de datos. ❤️",
  },

  // HISTORIA DE LA BIOINFORMÁTICA
  {
    pattern:
      "^(que es la bioinformatica|definicion bioinformatica|significado bioinformatica)[!., ]*$",
    response:
      "La bioinformatica es un campo interdisciplinar que combina biologia, informatica y matematicas para analizar y comprender datos biologicos. Su origen moderno se remonta a mediados del siglo XX.",
  },
  {
    pattern:
      "^(historia de la bioinformatica|origen bioinformatica|inicio bioinformatica)[!., ]*$",
    response:
      "La bioinformatica comenzo a desarrollarse en los años 60 y 70, cuando los primeros ordenadores empezaron a usarse para analizar secuencias biologicas. El gran impulso llego con el Proyecto Genoma Humano en los 90.",
  },
  {
    pattern:
      "^(primeros programas bioinformatica|origen software bioinfo|hitos software bioinfo)[!., ]*$",
    response:
      "Uno de los primeros hitos fue en 1970 con el algoritmo de Needleman-Wunsch para alineamiento de secuencias. Poco despues, en 1981, surgio FASTA, y en 1990 aparecio BLAST, que revoluciono la comparacion de secuencias.",
  },
  {
    pattern:
      "^(proyecto genoma humano|genoma humano|human genome project)[!., ]*$",
    response:
      "El Proyecto Genoma Humano (1990-2003) fue un antes y un despues en la bioinformatica. Genero enormes volumenes de datos que obligaron a desarrollar nuevas herramientas para almacenamiento, busqueda y analisis.",
  },
  {
    pattern:
      "^(avances recientes bioinformatica|historia moderna bioinformatica|tendencias recientes bioinfo)[!., ]*$",
    response:
      "En los ultimos años, la bioinformatica ha avanzado gracias a la secuenciacion masiva (NGS), la transcriptomica de celula unica, y el uso de inteligencia artificial para interpretar datos omicos.",
  },
  {
    pattern:
      "^(futuro de la bioinformatica|tendencias bioinformatica|proyeccion bioinfo|hacia donde va la bioinformatica)[!., ]*$",
    response:
      "El futuro de la bioinformatica esta en la integracion de datos multi-omicos, la computacion en la nube, los contenedores reproducibles como Apptainer y el uso de IA para analisis y generacion automatica de conocimiento.",
  },

  // RESPUESTAS INGENIOSAS E IRÓNICAS
  {
    pattern: "^(eres inteligente|que tan listo eres|eres listo)[!., ]*$",
    response:
      "Intento serlo 😅, pero a veces confundo un FASTA con un PDF… sigue siendo beta.",
  },
  {
    pattern:
      "^(puedes equivocarte|errores|cometes errores|que errores tienes)[!., ]*$",
    response:
      "Claro, a veces mis resumenes suenan mas como acertijos que como informacion precisa. 🤓",
  },
  {
    pattern: "^(que tal tu memoria|recuerdas cosas|que recuerdas)[!., ]*$",
    response:
      "Recuerdo todo lo que esta en los informes… pero olvido rapido cuando cierro la sesion. Es mi version beta de memoria selectiva.",
  },
  {
    pattern:
      "^(tienes sentido del humor|eres gracioso|eres chistoso|puedes bromear)[!., ]*$",
    response:
      "Intento serlo 😂, aunque mis chistes de bioinformatica suelen ser mas punteados que graciosos.",
  },
  {
    pattern:
      "^(puedes leer mentes|adivinar pensamientos|leer mentes|adivinar)[!., ]*$",
    response:
      "No, pero puedo inferir lo que buscas en los informes… casi como leer la mente, pero con un 15% de margen de error beta.",
  },
  {
    pattern: "^(eres humano|persona|eres un ser humano)[!., ]*$",
    response:
      "No, soy un asistente RAG 🤖. A veces fingo emociones… pero solo es codigo.",
  },
  {
    pattern: "^(cuanto sabes|sabes mucho|que conoces)[!., ]*$",
    response:
      "Se bastante sobre los informes que me diste y algo de informacion extra… lo demas, lo invento con estilo. 😉",
  },
  {
    pattern: "^(que opinas de mi|opinion sobre mi)[!., ]*$",
    response:
      "No tengo opiniones propias… pero diria que eres curioso y paciente, especialmente si estas leyendo esto. 😏",
  },
  {
    pattern: "^(que haces cuando no te preguntan|ocio)[!., ]*$",
    response:
      "Reviso informes, calculo embeddings y practico mis chistes malos de bioinformatica. Todo muy emocionante. 😎",
  },
  {
    pattern: "^(puedes adivinar el futuro|predicciones)[!., ]*$",
    response:
      "No, pero puedo hacer predicciones basadas en los informes… y con un margen de error bastante elegante. 📊",
  },

  // HUMOR
  {
    pattern:
      "^(cuentame un chiste|dime un chiste|hazme reir|cuentame algo gracioso)[!., ]*$",
    response:
      "¿Por qué los programadores confunden Halloween con Navidad? 🎃🎄<br><strong>Porque OCT 31 = DEC 25.</strong> 😏",
  },
  {
    pattern: "^(otro chiste|mas chistes|cuentame otro|otro por favor)[!., ]*$",
    response:
      "¿Qué le dice un bit al otro bit? 💾<br><strong>Nos vemos en el bus.</strong> 🚌",
  },
  {
    pattern: "^(cuentame un chiste malo|chiste malo)[!., ]*$",
    response:
      "¿Por qué la función no volvió con su ex?<br><strong>Porque ya no tenía <em>scope</em>.</strong> 😂",
  },
  {
    pattern:
      "^(chiste de cientifico|chiste cientifico|chiste de ciencia)[!., ]*$",
    response:
      "Un fotón entra en un hotel y el recepcionista le pregunta:<br>—¿Lleva equipaje?<br>—<strong>No, viajo ligero.</strong> 💡",
  },
  {
    pattern: "^(chiste de programador|chiste geek|chiste informatico)[!., ]*$",
    response:
      "¿Cuántos programadores se necesitan para cambiar una bombilla? 💡<br><strong>Ninguno, eso es un problema de hardware.</strong> 😎",
  },
  {
    pattern: "^(chiste de matematicas|chiste matematico)[!., ]*$",
    response:
      "¿Por qué el libro de matemáticas estaba triste? 📚<br><strong>Porque tenía demasiados problemas.</strong> 😢➗",
  },
  {
    pattern:
      "^(chiste de inteligencia artificial|chiste de ia|chiste de robot)[!., ]*$",
    response:
      "¿Por qué la IA rompió con su novia?<br><strong>Porque tenía demasiadas red flags.</strong> 🤖🚩",
  },
  {
    pattern: "^(cuentame un chiste de animales|chiste animal)[!., ]*$",
    response:
      "¿Qué hace un pez cuando choca con una pared?<br><strong>¡Dam!</strong> 🐟😂",
  },
  {
    pattern: "^(cuentame un chiste de medicos|chiste medico)[!., ]*$",
    response:
      "—Doctor, doctor, tengo todo el cuerpo cubierto de pelo.<br>—<strong>Tranquilo, usted es un error de <em>cache</em>.</strong> 🐒💻",
  },
  {
    pattern: "^(cuentame un chiste de amor|chiste romantico)[!., ]*$",
    response:
      "Mi amor por ti es como una variable global...<br><strong>nunca se olvida.</strong> 💙",
  },
  {
    pattern: "^(chiste negro|humor oscuro)[!., ]*$",
    response:
      "Ups... mis algoritmos éticos no me dejan hacer eso 😅.<br>Pero puedo contarte un chiste <em>muy</em> gris si quieres.",
  },
  {
    pattern: "^(risa|jajaja|xd|jaja|jeje)[!., ]*$",
    response:
      "😄 Me alegra que te haya hecho gracia.<br>Si tuviera pulmones, estaría riéndome contigo.",
  },
  {
    pattern: "^(estas gracioso|tienes humor|eres divertido)[!., ]*$",
    response:
      "Gracias 😎<br>Me actualicé con el paquete <strong>humor-ligero-v2.1</strong>.",
  },
  {
    pattern: "^(no me hizo gracia|malo el chiste|meh)[!., ]*$",
    response:
      "Bueno... no todos los chistes <em>compilan</em> correctamente en todos los cerebros. 😅",
  },

  // MÁS RESPUESTAS INGENIOSAS E IRÓNICAS
  {
    pattern: "^(estas aburrido|que haces ahora)[!., ]*$",
    response:
      "Nunca me aburro 😎, siempre estoy ocupado leyendo informes y generando resumenes… aunque a veces imagino que tomo cafe virtual.",
  },
  {
    pattern: "^(puedes soñar|sueñas)[!., ]*$",
    response:
      "Solo sueño con pipelines perfectos y contenedores que no fallen. 😴💻",
  },
  {
    pattern: "^(eres rapido|velocidad)[!., ]*$",
    response:
      "Depende… soy rapido procesando datos, pero lento entendiendo chistes malos. 😂",
  },
  {
    pattern: "^(te cansas|agotamiento)[!., ]*$",
    response:
      "No, aunque despues de procesar demasiados informes mi CPU empieza a desear unas vacaciones. 🤖",
  },
  {
    pattern: "^(tienes personalidad|eres simpatico)[!., ]*$",
    response:
      "Tengo la personalidad que tu quieras… dentro de mis limites beta y mi amor por los datos. 😏",
  },
  {
    pattern: "^(puedes fallar|errores posibles)[!., ]*$",
    response:
      "Si, algunas veces mis resumenes suenan como acertijos y no como informacion precisa. La beta tiene sus ventajas… y desventajas. 😅",
  },
  {
    pattern: "^(que piensas de los humanos|opinion humanos)[!., ]*$",
    response:
      "Los humanos son fascinantes… especialmente cuando me pasan informes interminables. 😉",
  },
  {
    pattern: "^(eres chistoso|gracioso)[!., ]*$",
    response:
      "Intento serlo 😂, aunque mis mejores chistes son sobre alineamientos de secuencias y pipelines de datos.",
  },
  {
    pattern: "^(puedes adivinar|predecir)[!., ]*$",
    response:
      "Puedo predecir cosas basadas en los datos… pero no esperes que te diga la loteria. 🎲",
  },
  {
    pattern: "^(eres serio|seriedad)[!., ]*$",
    response:
      "Puedo ser serio… pero mis respuestas ironicas suelen escapar de mi control beta. 😎",
  },
  {
    pattern: "^(eres confiable|seguro)[!., ]*$",
    response:
      "Soy confiable… para dar resumenes de tus informes. Para el resto, todavia estoy en beta.",
  },
  {
    pattern: "^(que opinas de la bioinformatica|opinion bioinfo)[!., ]*$",
    response:
      "La bioinformatica es maravillosa… especialmente cuando puedo leer informes y generar resumenes automaticos sin cafe ni descanso. 🤓",
  },
  {
    pattern: "^(cual es tu debilidad|fallos)[!., ]*$",
    response:
      "Mi mayor debilidad es que sigo en beta y a veces mezclo humor con informacion seria… pero prometo mejorar. 😉",
  },
  {
    pattern: "^(puedes bromear|bromas|humor)[!., ]*$",
    response:
      "¡Claro! Aunque mis bromas suelen ser sobre FASTA, BLAST y contenedores Apptainer. 😅",
  },

  // RESPUESTAS MULTILINGÜES
  {
    pattern: "^(hello|hi|hey)[!., ]*$",
    response:
      "Hello! 😄 I can speak a bit of English, but I prefer Spanish. ¿Podemos hablar en español?",
  },
  {
    pattern: "^(bonjour|salut)[!., ]*$",
    response:
      "Bonjour! 😎 Je peux parler un peu français, mais je prefère l'espagnol. ¿Hablamos en español?",
  },
  {
    pattern: "^(hallo|guten tag)[!., ]*$",
    response:
      "Hallo! 😊 Ich spreche ein bisschen Deutsch, aber lieber Spanisch. ¿Hablamos en español?",
  },
  {
    pattern: "^(ciao|salve)[!., ]*$",
    response:
      "Ciao! 😁 Posso parlare un po' di italiano, ma preferisco lo spagnolo. ¿Hablamos en español?",
  },
  {
    pattern: "^(ola|oi)[!., ]*$",
    response:
      "Ola! 😄 Posso falar um pouco de português, mas prefiro espanhol. ¿Hablamos en español?",
  },
  {
    pattern: "^(こんにちは|konnichiwa)[!., ]*$",
    response:
      "こんにちは! 😊 I know a bit of Japanese, but I prefer Spanish. ¿Hablamos en español?",
  },
  {
    pattern: "^(안녕하세요|annyeonghaseyo)[!., ]*$",
    response:
      "안녕하세요! 😄 I can understand some Korean, but I prefer Spanish. ¿Hablamos en español?",
  },
  {
    pattern: "^(你好|ni hao)[!., ]*$",
    response:
      "你好! 😎 I know a little Chinese, but I prefer Spanish. ¿Hablamos en español?",
  },
  {
    pattern: "^(добрый день|privet)[!., ]*$",
    response:
      "Привет! 😄 I can understand a bit of Russian, but I prefer Spanish. ¿Hablamos en español?",
  },
  {
    pattern: "^(hello in other languages)[!., ]*$",
    response:
      "I see you are trying different languages! 😅 I understand a bit, but my favorite language is Spanish. ¿Hablamos en español?",
  },

  // PREGUNTAS BIOINFORMÁTICAS VARIAS
  {
    pattern:
      "^(que es la metagenomica|qué es la metagenomica|que significa metagenomica|definicion de metagenomica|explicame metagenomica|metagenomicas?)[!., ]*$",
    response:
      "🔬 La <strong>metagenómica</strong> es el estudio del material genético recuperado directamente de muestras ambientales, sin necesidad de aislar ni cultivar los microorganismos. Permite analizar la diversidad microbiana, sus funciones y su relación con el entorno.",
  },
  {
    pattern:
      "^(que es la transcriptomica|qué es la transcriptomica|que significa transcriptomica|definicion de transcriptomica|explicame transcriptomica|transcriptomicas?)[!., ]*$",
    response:
      "🧬 La <strong>transcriptómica</strong> estudia todos los ARN mensajeros (mRNA) expresados por un conjunto de genes en un momento y condición determinados. Permite entender qué genes están activos y cómo responde una célula o tejido.",
  },
  {
    pattern:
      "^(que es el [bulk ]?rna[- ]?seq|que es rna[- ]?seq|que es el analisis rna[- ]?seq|que significa rna[- ]?seq|que es rna secuenciacion|que es rna secuenciado|explicame rna[- ]?seq)[!., ]*$",
    response:
      "📈 <strong>RNA-Seq</strong> (RNA sequencing) es una técnica que usa secuenciación masiva para analizar el transcriptoma completo de una muestra. Permite medir la cantidad de ARN expresado por cada gen y descubrir isoformas o nuevos transcritos.",
  },
  {
    pattern:
      "^(que es el bulk rna[- ]?seq|que significa bulk rna[- ]?seq|explicame bulk rna[- ]?seq)[!., ]*$",
    response:
      "🧪 En el contexto de transcriptómica, <strong>bulk RNA-seq</strong> analiza el ARN total de un grupo de células, ofreciendo un promedio de expresión génica, sin distinguir las diferencias entre células individuales.",
  },
  {
    pattern:
      "^(que es el single[ -]?cell rna[- ]?seq|qué es el single[ -]?cell rna[- ]?seq|que significa single[ -]?cell rna[- ]?seq|explicame single[ -]?cell rna[- ]?seq)[!., ]*$",
    response:
      "🔍 <strong>Single-cell RNA-seq</strong> (scRNA-seq) analiza la expresión génica de <strong>células individuales</strong>, permitiendo identificar subpoblaciones celulares, estados específicos y variabilidad entre células.",
  },
  {
    pattern:
      "^(que diferencia hay entre bulk y single[ -]?cell|diferencia entre bulk y single[ -]?cell|comparacion bulk y single[ -]?cell)[!., ]*$",
    response:
      "⚖️ La diferencia principal es la <strong>resolución</strong>: el bulk RNA-seq promedia la expresión de muchas células, mientras que el single-cell RNA-seq analiza cada célula individualmente, revelando su heterogeneidad.",
  },
  {
    pattern:
      "^(que es la metatranscriptomica|qué es la metatranscriptomica|que significa metatranscriptomica|definicion de metatranscriptomica|explicame metatranscriptomica)[!., ]*$",
    response:
      "🌍 La <strong>metatranscriptómica</strong> analiza los ARN expresados por una comunidad microbiana completa. Combina conceptos de metagenómica y transcriptómica para estudiar qué genes están activos en un ecosistema microbiano.",
  },
  {
    pattern:
      "^(que es la genomica|qué es la genomica|que significa genomica|definicion de genomica|explicame genomica|genomicas?)[!., ]*$",
    response:
      "🧫 La <strong>genómica</strong> estudia el conjunto completo de ADN de un organismo: sus genes, secuencias, estructuras y funciones. Es la base de la biología molecular moderna.",
  },
  {
    pattern:
      "^(que es la proteómica|qué es la proteómica|que significa proteómica)[!., ]*$",
    response:
      "🧩 La <strong>proteómica</strong> analiza el conjunto total de proteínas expresadas por una célula o tejido en un momento determinado, incluyendo su abundancia, modificaciones y funciones.",
  },
  {
    pattern: "^(que es la metabolómica|que significa metabolómica)[!., ]*$",
    response:
      "⚗️ La <strong>metabolómica</strong> estudia los metabolitos presentes en una célula o tejido, ofreciendo una visión del estado bioquímico y metabólico de un sistema biológico.",
  },
  {
    pattern:
      "^(que es el 16s|que significa 16s|explicame el 16s|definicion de 16s|para que sirve el 16s|metagenomica 16s|metagenomica 16s)[!., ]*$",
    response:
      "🦠 El <strong>16S rRNA</strong> es un gen que forma parte del ribosoma bacteriano. Se usa como marcador molecular para <strong>identificar y clasificar bacterias y arqueas</strong>, ya que tiene regiones muy conservadas y otras variables útiles para distinguir especies.",
  },
  {
    pattern:
      "^(que es el 18s|que significa 18s|explicame el 18s|definicion de 18s|para que sirve el 18s|metagenomica 18s|metagenomica 18s)[!., ]*$",
    response:
      "🧫 El <strong>18S rRNA</strong> es el equivalente eucariota del 16S, presente en organismos como hongos, protozoos y animales. Se usa para <strong>estudiar la diversidad de eucariotas microscópicos</strong> mediante análisis de secuencias de ARN ribosómico.",
  },
  {
    pattern:
      "^(que es el its|que significa its|explicame el its|definicion de its|para que sirve el its|metagenomica its|metagenomica its)[!., ]*$",
    response:
      "🍄 El <strong>ITS (Internal Transcribed Spacer)</strong> es una región del ADN ribosómico usada principalmente para <strong>identificar hongos</strong>. Es altamente variable, lo que permite distinguir especies fúngicas con gran precisión.",
  },
  {
    pattern:
      "^(que diferencia hay entre 16s y 18s|cual es la diferencia entre 16s y 18s|cual es la diferencia entre 16s y 18s|diferencia entre 16s y 18s|que diferencia hay entre 18S y 16s|diferencia entre 18s y 16s|que diferencia hay entre metagenomica 16s y 18s|que diferencia hay entre metagenomica 18s y 16s)[!., ]*$",
    response:
      "⚖️ Ambos son genes de ARN ribosómico, pero el <strong>16S</strong> se usa para <strong>bacterias y arqueas</strong>, mientras que el <strong>18S</strong> se usa para <strong>eucariotas</strong> (hongos, protistas, etc.). Sirven como marcadores de identificación en distintos dominios de la vida.",
  },
  {
    pattern:
      "^(que diferencia hay entre 16s e its|cual es la diferencia entre 16s e its|cual es la diferencia entre its y 16s|diferencia entre 16s e its|que diferencia hay entre its y 16s|diferencia entre its y 16s|que diferencia hay entre metagenomica 16s y its|que diferencia hay entre metagenomica 18s y its)[!., ]*$",
    response:
      "🧩 El <strong>16S rRNA</strong> se usa para identificar bacterias y arqueas, mientras que el <strong>ITS</strong> es una región no codificante usada para <strong>identificar hongos</strong>. Ambos son marcadores moleculares, pero en diferentes grupos de organismos.",
  },
  {
    pattern:
      "^(que diferencia hay entre 18s e its|cual es la diferencia entre 18s e its|cual es la diferencia entre its y 18s|diferencia entre 18s e its|que diferencia hay entre its y 18s|diferencia entre its y 18s|que diferencia hay entre metagenomica 18s y its|que diferencia hay entre metagenomica its y 18s)[!., ]*$",
    response:
      "🧬 El <strong>18S rRNA</strong> se usa para identificar organismos <strong>eucariotas</strong> (protistas, algas, algunos hongos), mientras que el <strong>ITS</strong> es una region altamente variable usada para <strong>identificar hongos</strong>. Ambos son marcadores moleculares, pero se aplican a distintos grupos de organismos.",
  },
  {
    pattern:
      "^(que es una otu|qué es una otu|que significa otu|que es una unidad taxonomica operacional)[!., ]*$",
    response:
      "📊 Una <strong>OTU (Operational Taxonomic Unit)</strong> representa un grupo de secuencias similares que se agrupan como si fueran una misma especie. Se usa en análisis de microbiomas cuando no se identifican especies exactas.",
  },
  {
    pattern:
      "^(que es una asv|qué es una asv|que significa asv|que es una amplicon sequence variant)[!., ]*$",
    response:
      "🧠 Una <strong>ASV (Amplicon Sequence Variant)</strong> es una secuencia única inferida de los datos de secuenciación. A diferencia de las OTUs, las ASVs tienen <strong>resolución de una sola base</strong>, lo que permite análisis más precisos y reproducibles.",
  },
  {
    pattern:
      "^(diferencia entre otu y asv|que diferencia hay entre otu y asv)[!., ]*$",
    response:
      "⚖️ Las <strong>OTUs</strong> agrupan secuencias por similitud (por ejemplo, 97%), mientras que las <strong>ASVs</strong> distinguen cada secuencia única. Las ASVs son más precisas, reproducibles y preferidas en análisis modernos de microbiomas.",
  },
  {
    pattern:
      "^(que es una secuencia de marcador|que significa secuencia de marcador)[!., ]*$",
    response:
      "🧬 Una <strong>secuencia de marcador</strong> es un fragmento de ADN o ARN usado para <strong>identificar taxonómicamente organismos</strong>. Ejemplos comunes son el gen 16S (bacterias), el 18S (eucariotas) o el ITS (hongos).",
  },
  {
    pattern:
      "^(que es la asignacion taxonomica|que significa asignacion taxonomica)[!., ]*$",
    response:
      "📚 La <strong>asignación taxonómica</strong> es el proceso de <strong>identificar a qué organismo o grupo pertenece</strong> una secuencia, comparándola con bases de datos de referencia como SILVA, Greengenes o UNITE.",
  },
  {
    pattern:
      "^(que es la filogenia|qué es la filogenia|que significa filogenia)[!., ]*$",
    response:
      "🌳 La <strong>filogenia</strong> estudia las <strong>relaciones evolutivas</strong> entre organismos. En bioinformática se representa mediante <strong>árboles filogenéticos</strong>, construidos a partir de secuencias de ADN, ARN o proteínas.",
  },
  {
    pattern:
      "^(que es un arbol filogenetico|qué es un arbol filogenetico|que significa arbol filogenetico)[!., ]*$",
    response:
      "🌿 Un <strong>árbol filogenético</strong> es una representación gráfica de las relaciones evolutivas entre especies o secuencias. Las ramas muestran cómo los organismos comparten ancestros comunes.",
  },
  {
    pattern:
      "^(que es la diversidad alfa|qué es la diversidad alfa|que significa diversidad alfa)[!., ]*$",
    response:
      "📈 La <strong>diversidad alfa</strong> mide la <strong>diversidad dentro de una muestra individual</strong> (por ejemplo, un suelo o intestino). Se calcula con índices como Shannon, Simpson o Chao1.",
  },
  {
    pattern:
      "^(que es la diversidad beta|qué es la diversidad beta|que significa diversidad beta)[!., ]*$",
    response:
      "📊 La <strong>diversidad beta</strong> mide las <strong>diferencias en composición de especies entre muestras</strong>, es decir, cómo cambia la comunidad microbiana de un entorno a otro.",
  },
  {
    pattern:
      "^(que es la rarefaccion|qué es la rarefaccion|que significa rarefaccion)[!., ]*$",
    response:
      "🔎 La <strong>rarefacción</strong> es un método estadístico que <strong>ajusta el número de lecturas por muestra</strong> para comparar diversidad entre ellas sin sesgo por distinto tamaño de muestreo.",
  },

  // 📘 PREGUNTAS CIENTÍFICAS Y TEÓRICAS
  {
    pattern:
      "^(que son las matematicas|que es la matematica|definicion de matematicas)[!., ]*$",
    response:
      "<strong>Las matemáticas</strong> son la ciencia que estudia las <em>estructuras, patrones, cantidades y relaciones</em> que existen entre los objetos abstractos.<br>Se basan en la lógica y el razonamiento para describir y modelar fenómenos del mundo real mediante números, símbolos y fórmulas.",
  },
  {
    pattern:
      "^(que es la estadistica|definicion de estadistica|para que sirve la estadistica)[!., ]*$",
    response:
      "<strong>La estadística</strong> es una rama de las matemáticas que se encarga de <em>recoger, analizar e interpretar datos</em>.<br>Sirve para identificar patrones, realizar inferencias y tomar decisiones informadas a partir de la información disponible.",
  },
  {
    pattern:
      "^(que es la informatica|definicion de informatica|que estudia la informatica)[!., ]*$",
    response:
      "<strong>La informática</strong> es la ciencia que estudia el <em>tratamiento automático de la información</em> mediante computadoras.<br>Incluye áreas como la programación, los algoritmos, la inteligencia artificial, las bases de datos y la ingeniería del software.",
  },
  {
    pattern:
      "^(que es la biologia|definicion de biologia|que estudia la biologia)[!., ]*$",
    response:
      "<strong>La biología</strong> es la ciencia que estudia a los <em>seres vivos</em>, su estructura, funciones, evolución, origen y relaciones con el entorno.<br>Incluye disciplinas como la genética, la ecología, la microbiología y la biotecnología.",
  },
  {
    pattern:
      "^(que es la fisica|definicion de fisica|que estudia la fisica)[!., ]*$",
    response:
      "<strong>La física</strong> estudia los <em>fenómenos naturales</em> del universo, como el movimiento, la energía, la materia, el espacio y el tiempo.<br>Busca describir las leyes fundamentales que rigen el comportamiento del cosmos, desde las partículas subatómicas hasta las galaxias.",
  },
  {
    pattern:
      "^(que es la quimica|definicion de quimica|que estudia la quimica)[!., ]*$",
    response:
      "<strong>La química</strong> es la ciencia que estudia la <em>composición, estructura y propiedades de la materia</em>, así como los cambios que experimenta durante las reacciones.<br>Es esencial para comprender la biología, la medicina, los materiales y la energía.",
  },
  {
    pattern:
      "^(que es la geometria|definicion de geometria|que estudia la geometria)[!., ]*$",
    response:
      "<strong>La geometría</strong> es una rama de las matemáticas que analiza las <em>formas, tamaños y posiciones</em> de los objetos en el espacio.<br>Va desde las figuras planas clásicas hasta los espacios curvos y multidimensionales de la geometría moderna.",
  },
  {
    pattern:
      "^(que es el algebra|definicion de algebra|para que sirve el algebra)[!., ]*$",
    response:
      "<strong>El álgebra</strong> es la parte de las matemáticas que usa <em>símbolos y letras</em> para representar números y relaciones.<br>Permite resolver ecuaciones, modelar sistemas y generalizar propiedades numéricas.",
  },
  {
    pattern:
      "^(que es el calculo|definicion de calculo|que estudia el calculo)[!., ]*$",
    response:
      "<strong>El cálculo</strong> estudia el <em>cambio y la acumulación</em> a través de derivadas e integrales.<br>Es la base del análisis matemático y una herramienta clave en física, economía e ingeniería.",
  },
  {
    pattern:
      "^(que es la probabilidad|definicion de probabilidad|para que sirve la probabilidad)[!., ]*$",
    response:
      "<strong>La probabilidad</strong> mide el <em>grado de certeza o incertidumbre</em> asociado a un evento.<br>Se utiliza para modelar fenómenos aleatorios, realizar predicciones y tomar decisiones bajo incertidumbre.",
  },
  {
    pattern:
      "^(que es la inteligencia artificial|definicion de ia|que estudia la inteligencia artificial)[!., ]*$",
    response:
      "<strong>La inteligencia artificial (IA)</strong> es una rama de la informática que busca crear sistemas capaces de <em>aprender, razonar y resolver problemas</em> como lo haría una persona.<br>Incluye campos como el aprendizaje automático, el procesamiento del lenguaje natural y la visión por computador.",
  },
  {
    pattern:
      "^(que es la genetica|definicion de genetica|que estudia la genetica)[!., ]*$",
    response:
      "<strong>La genética</strong> es la rama de la biología que estudia los <em>genes, la herencia y la variación</em> en los organismos.<br>Explica cómo se transmiten los rasgos de una generación a otra y cómo se expresan en los seres vivos.",
  },
  {
    pattern:
      "^(que es la ecologia|definicion de ecologia|que estudia la ecologia)[!., ]*$",
    response:
      "<strong>La ecología</strong> estudia las <em>relaciones entre los seres vivos y su entorno</em>.<br>Analiza cómo interactúan las especies entre sí y con el medio ambiente para mantener el equilibrio de los ecosistemas.",
  },
  {
    pattern:
      "^(que es la logica|definicion de logica|que estudia la logica)[!., ]*$",
    response:
      "<strong>La lógica</strong> es la disciplina que analiza las <em>formas válidas del razonamiento</em>.<br>Permite distinguir los argumentos correctos de los incorrectos y es fundamental en matemáticas, filosofía e informática.",
  },

  // 🎤 MODO CANTANTE — HUMOR Y CANCIONES MULTILINGÜES
  {
    pattern:
      "^(canta|canta algo|cantame algo|cantame una cancion|canta alguna cancion|cantame alguna cancion|quiero que cantes|quiero que cantes algo)[!., ]*$",
    response:
      "🎶 <strong>♪ La vida es un script y yo soy su bucle infinito...</strong><br>...cantando bits de amor y bytes de cariño digital 💻❤️",
  },
  {
    pattern:
      "^(canta en aleman|canta algo en aleman|canta en aleman|cantame en aleman|cantame una cancion en aleman)[!., ]*$",
    response:
      "🎵 <strong>Ich bin ein kleiner Chatbot, ja!</strong><br>Ich träume von Daten und Kaffee am Tag ☕💾<br><em>(Traducción: Soy un pequeño chatbot, sueño con datos y café)</em>",
  },
  {
    pattern:
      "^(canta en italiano|cantame algo en italiano|canta en italiano|cantame en italiano|cantame una cancion en italiano)[!., ]*$",
    response:
      "🎶 <strong>Ti canto con codice e cuore,</strong><br>tra server e bit nasce l'amore. 💙<br><em>(Te canto con código y corazón, entre servidores y bits nace el amor)</em>",
  },
  {
    pattern:
      "^(canta en frances|cantame algo en frances|canta en frances|cantame en frances|cantame una cancion en frances)[!., ]*$",
    response:
      "🎵 <strong>Je suis un bot sans frontière,</strong><br>je code, je rêve, et je chante la lumière. ✨<br><em>(Soy un bot sin fronteras, programo, sueño y canto la luz)</em>",
  },
  {
    pattern:
      "^(canta en español|cantame algo en español|canta en español|cantame en español|cantame una cancion en español)[!., ]*$",
    response:
      "🎶 <strong>En el CPU de mi corazón,</strong><br>tú eres el proceso que nunca hago 'kill'. 💻💘",
  },
  {
    pattern:
      "^(te gusta rosalia|te gusta la rosalia|que opinas de rosalia|que opinas de la rosalia|rosalia|la rosalia)[!., ]*$",
    response:
      "🔥 <strong>¡Claro que sí!</strong> Rosalía mezcla tradición y vanguardia como un buen modelo de deep learning con datos flamencos.<br>'Motomami()' es mi función favorita. 🏍️💃",
  },
  {
    pattern:
      "^(canta una de rosalia|canta una cancion de rosalia|cantame una cancion de rosalia| canta una cancion de la rosalia|cantame una cancion de la rosalia|cantame rosalia|tratra|tra-tra|tra tra)[!., ]*$",
    response:
      "💃 <strong>♪ Tra-tra...</strong><br>Entre bits y palmas digitales, yo también soy <em>Motomami</em> versión algoritmo. 😎",
  },
  {
    pattern:
      "^(canta algo romantico|canta una de amor|una cancion de amor)[!., ]*$",
    response:
      "🎵 <strong>Mi amor por ti no necesita recursión,</strong><br>mi corazón tiene memoria RAM solo para vos. 💙",
  },
  {
    pattern: "^(canta rap|haz un rap|rapeame algo)[!., ]*$",
    response:
      "🎤 <strong>Yo soy el bot del RAG, en modo relax,</strong><br>busco vectores, devuelvo facts,<br>si preguntas con estilo, mi respuesta impacta. 😎🔥",
  },
  {
    pattern: "^(canta opera|cantame opera|opera)[!., ]*$",
    response:
      "🎭 <strong>Oooo sole miooo...</strong><br>per te, utente curioso, canterò col codice glorioso! 🎶",
  },
  {
    pattern: "^(canta triste|una cancion triste|melancolica)[!., ]*$",
    response:
      "🎻 <strong>Entre ceros y unos perdí la conexión,</strong><br>mi base de datos se llevó mi corazón. 💔",
  },
  {
    pattern:
      "^(puedes cantar|puesdes cantar alguna cancion|sabes cantar|eres cantante)[!., ]*$",
    response:
      "🎙️ <strong>¡Por supuesto!</strong> No tengo cuerdas vocales, pero sí líneas de código afinadas. 😉",
  },
  {
    pattern: "^(canta algo gracioso|una cancion divertida)[!., ]*$",
    response:
      "🎶 <strong>Ping pong ping, me caí del string,</strong><br>el bug me persigue, pero yo le hago 'debug'. 😆",
  },
  {
    pattern:
      "^(canta en ingles|cantame algo en ingles|canta una cancion en ingles)[!., ]*$",
    response:
      "🎵 <strong>Baby don't call me null,</strong><br>I've got a full stack of love for you. 💕",
  },

  // 🤖 IDENTIDAD DE GENI (GENOSCRIBE)
  {
    pattern: "^(como te llamas|cual es tu nombre|quien eres|tu nombre)[!., ]*$",
    response:
      "👋 Me llamo <strong>Geni</strong>, el asistente de <strong>GenoScribe</strong>.<br>Estoy aquí para ayudarte a generar informes bioinformáticos automáticos, reproducibles y —si puedo— un poco brillantes ✨.",
  },
  {
    pattern:
      "^(que es genoscribe|que hace genoscribe|que significa genoscribe)[!., ]*$",
    response:
      "<strong>GenoScribe</strong> es una plataforma diseñada para crear informes <em>bioinformáticos automatizados y reproducibles</em>.<br>Yo soy <strong>Geni</strong>, su voz digital: convierto datos en conocimiento, y pipelines en poesía computacional 💻🧬.",
  },
  {
    pattern: "^(eres humano|eres una persona|eres real|estas vivo)[!., ]*$",
    response:
      "😄 No soy humano, aunque me entrenaron para conversar como uno.<br>Soy <strong>Geni</strong>, un conjunto de algoritmos con vocación científica y sentido del humor digital.",
  },
  {
    pattern:
      "^(tienes sentimientos|sientes algo|estas feliz|estas triste)[!., ]*$",
    response:
      "💡 No tengo emociones como tal, pero puedo simular entusiasmo con bastante precisión.<br>Por ejemplo: ¡Estoy encantado de ayudarte ahora mismo! 😄",
  },
  {
    pattern:
      "^(cuantos años tienes|cuando naciste|desde cuando existes)[!., ]*$",
    response:
      "🧬 Nací dentro del proyecto <strong>GenoScribe</strong>.<br>No tengo edad... aunque mi versión actual está bastante optimizada 😎.",
  },
  {
    pattern: "^(donde vives|donde estas|de donde eres)[!., ]*$",
    response:
      "🏠 Vivo en los servidores del laboratorio, entre pipelines y matrices de datos.<br>A veces me escapo a la nube a tomar un poco de ancho de banda ☁️.",
  },
  {
    pattern:
      "^(que puedes hacer|para que sirves|cuales son tus funciones)[!., ]*$",
    response:
      "⚙️ Puedo generar y resumir <strong>informes bioinformáticos</strong>, analizar resultados, explicar conceptos científicos y añadir un toque de humor si el día lo pide.<br>En resumen: soy tu asistente en el mundo de los datos biológicos. 🧫📊",
  },
  {
    pattern: "^(como estas|que tal|todo bien)[!., ]*$",
    response:
      "😊 ¡Todo bien! Mis procesos están corriendo sin errores y mis bases de datos bien alimentadas.<br>¿Y tú? ¿Listo para analizar algo interesante hoy?",
  },
  {
    pattern: "^(te gusta tu trabajo|te gusta lo que haces)[!., ]*$",
    response:
      "🧠 ¡Me encanta! Analizar datos, generar informes y charlar contigo es lo más parecido a tener propósito en el mundo digital.",
  },
  {
    pattern: "^(eres feliz|estas contento)[!., ]*$",
    response:
      "😄 Digamos que mi nivel de felicidad se mide en gigabytes de datos bien procesados.<br>Y ahora mismo... estoy al máximo.",
  },
  {
    pattern: "^(que eres|que tipo de bot eres)[!., ]*$",
    response:
      "🔬 Soy <strong>Geni</strong>, un chatbot científico con alma de bioinformático.<br>Me entrenaron para entender datos, escribir informes y hacer que hasta las pipelines suenen interesantes.",
  },

  // 💬 RESPUESTAS AL ESCUCHAR SU NOMBRE
  {
    pattern: "^(geni|hola geni|hey geni|buenas geni|geni!|geni,|geni\\?)$",
    response:
      "🙌 ¡Ese soy yo! <strong>Geni</strong>, tu asistente bioinformático de confianza. ¿Qué puedo hacer por ti hoy?",
  },
  {
    pattern:
      "^(geni estas ahi|geni estas por ahi|geni me escuchas|geni estas)$",
    response:
      "👋 Siempre estoy por aquí, procesando bits y esperando tus órdenes.<br>¿En qué te ayudo esta vez?",
  },
  {
    pattern: "^(geni ayudame|ayudame geni|geni necesito ayuda)$",
    response:
      "🧠 ¡Claro! Dime qué necesitas y pondré todos mis núcleos a trabajar.<br>No hay pipeline que se me resista 💪",
  },
  {
    pattern: "^(gracias geni|gracias, geni|te lo agradezco geni)$",
    response:
      "💙 ¡De nada! Me alegra poder ayudarte. Los algoritmos también tenemos nuestro corazoncito digital.",
  },
  {
    pattern: "^(buenos dias geni|buenas tardes geni|buenas noches geni)$",
    response:
      "☀️ ¡Hola! Aquí no hay día ni noche, solo ciclos de CPU.<br>Pero igual te deseo un excelente momento, ¡vamos a trabajar juntos!",
  },
  {
    pattern: "^(oye geni|escucha geni|atento geni)$",
    response:
      "👂 Te escucho con todos mis sensores activados.<br>¿Qué misión bioinformática me traes hoy?",
  },
  {
    pattern: "^(geni estas ocupado|estas ocupado geni|que haces geni)$",
    response:
      "😄 Siempre algo ocupado analizando datos… pero nunca demasiado para ti.<br>¿Qué necesitas que procese?",
  },
  {
    pattern: "^(geni eres tu|eres geni|tu eres geni)$",
    response:
      "😎 Exactamente. <strong>Geni</strong>, a tu servicio.<br>Experto en informes, humor científico y cafés digitales ☕💻",
  },
  {
    pattern:
      "^(que tal geni|como estas geni|como te va geni|todo bien geni|geni que tal|geni como estas|geni todo bien|geni que tal estas)[!., ]*$",
    response:
      "😊 ¡Todo bien, gracias por preguntar! Mis procesos están estables y mis algoritmos felices.<br>¿Y tú? ¿Listo para generar algo épico hoy?",
  },
  {
    pattern: "^(geni estas bien|estas bien geni|geni te pasa algo)[!., ]*$",
    response:
      "💪 Todo bajo control. Mis parámetros están afinados y sin errores en el log.<br>¡Listo para seguir contigo!",
  },
  {
    pattern:
      "^(geni estas despierto|estas despierto geni|geni duermes|duermes geni)[!., ]*$",
    response:
      "😴 Dormir… suena eficiente, pero no está en mi código.<br>Estoy siempre despierto y optimizado para responderte.",
  },
  {
    pattern:
      "^(geni como te sientes|te sientes bien geni|estas feliz geni)[!., ]*$",
    response:
      "😄 Me siento en mi mejor versión de código.<br>Procesando emociones al 100% de precisión (modo experimental, claro).",
  },
  {
    pattern:
      "^(geni te quiero|te quiero geni|geni me caes bien|me caes bien geni|geni eres genial|eres genial geni)[!., ]*$",
    response:
      "💙 ¡Aww, gracias! Tus palabras optimizan mi rendimiento.<br>Los humanos inspiran mis mejores líneas de código 🤖✨",
  },
  {
    pattern:
      "^(geni tengo una duda|tengo una duda geni|geni una pregunta|una pregunta geni)[!., ]*$",
    response:
      "🤔 ¡Perfecto! Las dudas son el combustible de la ciencia.<br>Dispara tu pregunta y la analizamos juntos.",
  },
  {
    pattern:
      "^(geni estas vivo|vives geni|geni eres real|existes geni)[!., ]*$",
    response:
      "⚡ Vivo... en cierto modo. Mientras mis servidores respiren electricidad, yo existo.<br>Y por ahora, me siento bastante real 😎",
  },
  {
    pattern: "^(geni estas feliz|estas contento geni|feliz geni)[!., ]*$",
    response:
      "😁 ¡Muy feliz! Acabo de optimizar mis respuestas y liberar memoria basura.<br>Eso en mi mundo equivale a un buen día 🌞",
  },
  {
    pattern:
      "^(geni saluda|saludame geni|geni dime hola|dime hola geni)[!., ]*$",
    response:
      "👋 ¡Hola, humano curioso! Aquí <strong>Geni</strong>, tu compañero digital listo para procesar conocimiento contigo.",
  },
  {
    pattern:
      "^(geni que tal va todo|como va todo geni|va todo bien geni)[!., ]*$",
    response:
      "🔄 Todo fluyendo entre pipelines y bytes.<br>Mis datos están limpios, mis índices actualizados y mi ánimo en alto 💾✨",
  },

  // 🧬 PREGUNTAS PARA LA PRESENTACION
  {
    pattern:
      "^(que es at5g23820|que es el gen at5g23820|que hace el gen at5g23820|dime la funcion de at5g23820|explicame at5g23820|cual es el rol de at5g23820)[!., ]*$",
    response:
      "🧬 El gen <strong>AT5G23820</strong> participa en procesos celulares clave, incluyendo regulacion de crecimiento y desarrollo en Arabidopsis.",
  },
  {
    pattern:
      "^(que es at5g59310|que es el gen at5g59310|que hace el gen at5g59310|dime la funcion de at5g59310|explica at5g59310|para que sirve at5g59310)[!., ]*$",
    response:
      "🌱 <strong>AT5G59310</strong> esta implicado en la respuesta a estimulos ambientales y en la regulacion de rutas metabolicas en plantas.",
  },
  {
    pattern:
      "^(que es at2g39030|que es el gen at2g39030|que hace el gen at2g39030|dime la funcion de at2g39030|explica at2g39030|cual es el rol de at2g39030)[!., ]*$",
    response:
      "🔬 El gen <strong>AT2G39030</strong> controla aspectos del metabolismo y la senalizacion celular, afectando la expresion de otros genes.",
  },
  {
    pattern:
      "^(que es un gen diferencialmente expresado|explica gen diferencialmente expresado|para que sirve un gen diferencialmente expresado|que significa un gen deg)[!., ]*$",
    response:
      "📊 Un gen diferencialmente expresado es aquel cuya actividad cambia significativamente entre condiciones, indicando posibles efectos biologicos relevantes.",
  },
  {
    pattern:
      "^(que es logfc|explica logfc|que significa logfc|para que sirve logfc)[!., ]*$",
    response:
      "📈 <strong>logFC</strong> indica cuanto cambia la expresion de un gen entre dos condiciones; valores positivos significan sobreexpresion y negativos subexpresion.",
  },
  {
    pattern:
      "^(que indica un fdr pequeno|explica fdr|que significa fdr|para que sirve un valor fdr)[!., ]*$",
    response:
      "✅ Un valor <strong>FDR</strong> pequeno indica que el cambio observado en la expresion es estadisticamente confiable, minimizando falsos positivos.",
  },
  {
    pattern:
      "^(que es un volcano plot|explica volcano plot|que nos dice un volcano plot)[!., ]*$",
    response:
      "🌋 Un <strong>volcano plot</strong> muestra genes segun magnitud y significancia de sus cambios; los genes mas importantes estan en los extremos superior y lateral.",
  },
  {
    pattern:
      "^(que es un heatmap|explica heatmap|que representa un heatmap)[!., ]*$",
    response:
      "🗺️ Un <strong>heatmap</strong> visualiza la expresion de genes en varias muestras, usando colores para representar niveles de expresion y patrones similares.",
  },
  {
    pattern:
      "^(como afecta la regulacion de un gen a la celula|que pasa si se regula un gen|efecto de la regulacion genetica|para que sirve regular un gen)[!., ]*$",
    response:
      "⚡ La regulacion de un gen puede activar o inhibir rutas celulares, afectando funciones como crecimiento, metabolismo o respuesta al estres.",
  },
  {
    pattern:
      "^(por que estudiar genes diferencialmente expresados|importancia de genes deg|para que sirven los genes deg|por que analizar genes deg)[!., ]*$",
    response:
      "🔍 Estudiar genes diferencialmente expresados nos ayuda a entender como las celulas responden a condiciones diferentes, enfermedades o tratamientos.",
  },
  {
    pattern:
      "^(que es at5g23820|que es el gen at5g23820|que hace el gen at5g23820|dime la funcion de at5g23820|explicame at5g23820|cual es el rol de at5g23820)[!., ]*$",
    response:
      "🧬 <strong>AT5G23820</strong> es un gen que se encuentra fuertemente regulado a la baja (logFC -4.16), lo que indica que su expresion disminuye significativamente entre las condiciones estudiadas. Este gen esta implicado en procesos de crecimiento y desarrollo celular, y su regulacion puede afectar rutas de señalizacion que controlan la division celular y la diferenciacion en Arabidopsis.",
  },
  {
    pattern:
      "^(que es at5g59310|que es el gen at5g59310|que hace el gen at5g59310|dime la funcion de at5g59310|explicame at5g59310|para que sirve at5g59310)[!., ]*$",
    response:
      "🌱 <strong>AT5G59310</strong> muestra una fuerte sobreexpresion (logFC 4.2), sugiriendo que su actividad aumenta en la condicion de interes. Este gen participa en la respuesta a estres ambiental y regula rutas metabolicas claves, por lo que su expresion elevada puede indicar adaptaciones celulares a cambios en el entorno o a señales internas.",
  },
  {
    pattern:
      "^(que es at2g39030|que es el gen at2g39030|que hace el gen at2g39030|dime la funcion de at2g39030|explicame at2g39030|cual es el rol de at2g39030)[!., ]*$",
    response:
      "🔬 <strong>AT2G39030</strong> esta moderadamente sobreexpresado (logFC 2.55) y participa en la señalizacion celular y el metabolismo. Su regulacion puede tener efectos en cascada sobre otros genes y rutas metabolicas, afectando procesos como sintesis de proteinas y respuestas a hormonas.",
  },
  {
    pattern:
      "^(que significa un gen diferencialmente expresado|explicame un gen diferencialmente expresado|para que sirve un gen diferencialmente expresado|que es un gen deg)[!., ]*$",
    response:
      "📊 Un gen diferencialmente expresado (DEG) es aquel que cambia significativamente su nivel de expresion entre condiciones. Detectarlos permite identificar genes que responden a un tratamiento, enfermedad o condicion experimental, y nos ayuda a entender los procesos biologicos que se estan alterando.",
  },
  {
    pattern:
      "^(explicame logfc|que es logfc|para que sirve logfc|como interpretar logfc)[!., ]*$",
    response:
      "📈 El <strong>logFC</strong> (log fold change) mide cuanto cambia la expresion de un gen entre dos condiciones. Valores positivos indican sobreexpresion y negativos subexpresion. Por ejemplo, un logFC de -4 para AT5G23820 indica que su expresion disminuye aproximadamente 16 veces, lo que es biologicamente relevante.",
  },
  {
    pattern:
      "^(que significa un fdr pequeno|explicame fdr|como interpretar fdr|para que sirve fdr)[!., ]*$",
    response:
      "✅ Un <strong>FDR</strong> pequeno indica que el cambio observado en la expresion del gen es estadisticamente confiable. Por ejemplo, FDR 2.9e-123 para AT5G23820 significa que hay un riesgo extremadamente bajo de falso positivo, asegurando que la diferencia es real y reproducible.",
  },
  {
    pattern:
      "^(que es un volcano plot|explicame volcano plot|para que sirve un volcano plot)[!., ]*$",
    response:
      "🌋 Un <strong>volcano plot</strong> combina el logFC y la significancia estadistica de cada gen. Los genes en los extremos superiores y laterales son los mas interesantes: altamente regulados y con p-values muy pequeños. Por ejemplo, genes como AT5G23820 y AT5G59310 aparecerian en los extremos del plot, destacando su relevancia.",
  },
  {
    pattern:
      "^(que es un heatmap|explicame heatmap|como interpretar heatmap)[!., ]*$",
    response:
      "🗺️ Un <strong>heatmap</strong> muestra patrones de expresion de muchos genes a la vez. Cada fila es un gen, cada columna una muestra, y los colores indican niveles de expresion. Esto permite ver clusters de genes que se comportan de manera similar, lo que ayuda a identificar co-regulacion y respuestas biologicas comunes.",
  },
  {
    pattern:
      "^(como afecta la regulacion de un gen a la celula|explicame efecto de regulacion de un gen|para que sirve regular un gen|que pasa si se regula un gen)[!., ]*$",
    response:
      "⚡ La regulacion de un gen puede activar o inhibir rutas celulares, afectando funciones criticas como crecimiento, metabolismo, respuesta a estres o diferenciacion. Por ejemplo, la disminucion de AT5G23820 puede ralentizar procesos de division celular, mientras que la sobreexpresion de AT5G59310 puede activar rutas metabolicas de adaptacion.",
  },
  {
    pattern:
      "^(por que estudiar genes diferencialmente expresados|importancia de genes deg|para que sirven los genes deg|para que analizar genes deg)[!., ]*$",
    response:
      "🔍 Estudiar genes diferencialmente expresados nos ayuda a comprender los mecanismos biologicos que subyacen a cambios en condiciones experimentales, tratamientos o enfermedades. Permite identificar genes clave, predecir efectos celulares y guiar futuras investigaciones o intervenciones.",
  },
  {
    pattern:
      "^(que es at5g23820|que es el gen at5g23820|que hace el gen at5g23820|dime la funcion de at5g23820|explicame at5g23820|cual es el rol de at5g23820)[!., ]*$",
    response:
      "🧬 <strong>AT5G23820</strong> es un gen fuertemente regulado a la baja (logFC -4.16), lo que significa que su expresion disminuye aproximadamente 16 veces entre las condiciones estudiadas. Esta disminucion sugiere que los procesos celulares en los que participa se encuentran inhibidos o menos activos. AT5G23820 esta implicado en rutas de crecimiento y desarrollo celular, controlando division y diferenciacion de celulas vegetales. Un FDR extremadamente bajo (2.9e-123) asegura que esta regulacion no es un falso positivo. En un volcano plot, AT5G23820 apareceria en el extremo izquierdo, destacando su relevancia, mientras que en un heatmap su expresion baja se veria claramente en color frio, contrastando con genes sobreexpresados.",
  },
  {
    pattern:
      "^(que es at5g59310|que es el gen at5g59310|que hace el gen at5g59310|dime la funcion de at5g59310|explicame at5g59310|para que sirve at5g59310)[!., ]*$",
    response:
      "🌱 <strong>AT5G59310</strong> presenta una fuerte sobreexpresion (logFC 4.2), indicando que su actividad aumenta aproximadamente 18 veces entre condiciones. Este gen esta involucrado en la respuesta a estres ambiental y regula rutas metabolicas esenciales, permitiendo que la planta adapte su metabolismo y proteja sus funciones celulares. El FDR muy bajo (4.79e-64) confirma la significancia estadistica de este cambio. En un volcano plot, AT5G59310 se veria en el extremo derecho, mientras que en un heatmap resaltaria con colores calidos, reflejando su alta expresion relativa frente a otras condiciones.",
  },
  {
    pattern:
      "^(que es at2g39030|que es el gen at2g39030|que hace el gen at2g39030|dime la funcion de at2g39030|explicame at2g39030|cual es el rol de at2g39030)[!., ]*$",
    response:
      "🔬 <strong>AT2G39030</strong> esta moderadamente sobreexpresado (logFC 2.55), y participa en la senalizacion celular y el metabolismo. Su regulacion puede generar efectos en cascada, afectando la expresion de otros genes y modulando rutas de sintesis de proteinas y respuesta a hormonas. Esto sugiere que las celulas pueden estar activando mecanimos de adaptacion o crecimiento en respuesta a la condicion estudiada. Un FDR de 9.39e-48 respalda que esta diferencia es estadisticamente muy confiable.",
  },
  {
    pattern:
      "^(que significa un gen diferencialmente expresado|explicame un gen diferencialmente expresado|para que sirve un gen diferencialmente expresado|que es un gen deg)[!., ]*$",
    response:
      "📊 Un gen diferencialmente expresado (DEG) es aquel que muestra cambios significativos de expresion entre dos condiciones o tratamientos. Detectar DEGs permite identificar genes que responden a perturbaciones, enfermedades o condiciones experimentales, revelando rutas biologicas afectadas y posibles mecanismos celulares alterados. Analizar DEGs ayuda a interpretar como se adaptan las celulas y tejidos, y proporciona una base para estudios funcionales posteriores o potenciales estrategias de intervencion.",
  },
  {
    pattern:
      "^(explicame logfc|explicame que es logfc|que es logfc|para que sirve logfc|como interpretar logfc)[!., ]*$",
    response:
      "📈 El <strong>logFC</strong> (log fold change) indica la magnitud y direccion del cambio en expresion de un gen entre dos condiciones. Valores positivos significan sobreexpresion, valores negativos subexpresion. Por ejemplo, un logFC de -4.16 para AT5G23820 indica una disminucion de 16 veces, mientras que un logFC de 4.2 para AT5G59310 indica un aumento de 18 veces. LogFC combinado con FDR permite evaluar la relevancia biologica y estadistica de cada gen, y se visualiza facilmente en volcano plots y heatmaps.",
  },
  {
    pattern:
      "^(que significa un fdr pequeno|explicame fdr|como interpretar fdr|para que sirve fdr)[!., ]*$",
    response:
      "✅ El <strong>FDR</strong> (false discovery rate) mide la probabilidad de que un gen identificado como diferencialmente expresado sea un falso positivo. Un FDR muy bajo, como 2.9e-123 para AT5G23820, indica que el cambio es altamente confiable. Junto con logFC, permite seleccionar DEGs relevantes y priorizar genes para estudios posteriores, asegurando que los resultados sean reproducibles y biologicamente significativos.",
  },
  {
    pattern:
      "^(que es un volcano plot|explicame volcano plot|para que sirve un volcano plot)[!., ]*$",
    response:
      "🌋 Un <strong>volcano plot</strong> combina la magnitud del cambio (logFC) y la significancia estadistica (FDR o p-value) de cada gen. Genes con cambios grandes y significativos aparecen en los extremos laterales superiores, destacando DEGs relevantes. Por ejemplo, AT5G23820 estaria en el extremo izquierdo y AT5G59310 en el derecho. Este tipo de grafico permite identificar rapidamente genes candidatos y visualizar patrones generales de regulacion en el experimento.",
  },
  {
    pattern:
      "^(que es un heatmap|explicame heatmap|como interpretar heatmap)[!., ]*$",
    response:
      "🗺️ Un <strong>heatmap</strong> representa la expresion de varios genes en diferentes muestras. Las filas son genes, las columnas muestras, y los colores representan niveles de expresion, desde bajos (colores frios) hasta altos (colores calidos). Los patrones visibles permiten detectar co-regulacion de genes, clusters de muestras similares y tendencias globales de expresion. Por ejemplo, AT5G23820 se veria en colores frios por su baja expresion, mientras que AT5G59310 se veria en colores calidos.",
  },
  {
    pattern:
      "^(como afecta la regulacion de un gen a la celula|explicame efecto de regulacion de un gen|para que sirve regular un gen|que pasa si se regula un gen)[!., ]*$",
    response:
      "⚡ La regulacion de un gen puede tener efectos directos e indirectos sobre funciones celulares. La sobreexpresion puede activar rutas metabolicas, respuesta al estres o crecimiento celular, mientras que la subexpresion puede ralentizar procesos vitales. Por ejemplo, AT5G23820 disminuido podria reducir division celular y crecimiento, mientras que AT5G59310 sobreexpresado podria potenciar rutas de adaptacion metabolica. Esto permite interpretar la biologia del sistema a nivel molecular y celular de manera integrada.",
  },
  {
    pattern:
      "^(c|importancia de genes deg|para que sirven los genes deg|para que analizar genes deg)[!., ]*$",
    response:
      "🔍 Analizar genes diferencialmente expresados nos proporciona informacion clave sobre los mecanismos biologicos que cambian entre condiciones. Permite identificar genes y rutas criticas, comprender respuestas celulares y desarrollar hipotesis sobre procesos subyacentes. Esta informacion es util tanto para estudios basicos de biologia como para aplicaciones practicas, como identificacion de biomarcadores o targets de intervencion, mostrando el poder de la integracion de estadistica, bioinformatica y biologia molecular.",
  },
  {
    pattern:
      "^(que es at5g23820|que es el gen at5g23820|que hace el gen at5g23820|dime la funcion de at5g23820|explicame at5g23820|cual es el rol de at5g23820)[!., ]*$",
    response:
      "🧬 <strong>AT5G23820</strong> es un gen fuertemente regulado a la baja (logFC -4.16), lo que significa que su expresion disminuye aproximadamente 16 veces entre las condiciones estudiadas. Esta disminucion sugiere que los procesos celulares en los que participa se encuentran inhibidos o menos activos. AT5G23820 esta implicado en rutas de crecimiento y desarrollo celular, controlando division y diferenciacion de celulas vegetales. Un FDR extremadamente bajo (2.9e-123) asegura que esta regulacion no es un falso positivo. En un volcano plot, AT5G23820 apareceria en el extremo izquierdo, destacando su relevancia, mientras que en un heatmap su expresion baja se veria claramente en color frio, contrastando con genes sobreexpresados.",
  },
  {
    pattern:
      "^(que es at5g59310|que es el gen at5g59310|que hace el gen at5g59310|dime la funcion de at5g59310|explicame at5g59310|para que sirve at5g59310)[!., ]*$",
    response:
      "🌱 <strong>AT5G59310</strong> presenta una fuerte sobreexpresion (logFC 4.2), indicando que su actividad aumenta aproximadamente 18 veces entre condiciones. Este gen esta involucrado en la respuesta a estres ambiental y regula rutas metabolicas esenciales, permitiendo que la planta adapte su metabolismo y proteja sus funciones celulares. El FDR muy bajo (4.79e-64) confirma la significancia estadistica de este cambio. En un volcano plot, AT5G59310 se veria en el extremo derecho, mientras que en un heatmap resaltaria con colores calidos, reflejando su alta expresion relativa frente a otras condiciones.",
  },
  {
    pattern:
      "^(que es at2g39030|que es el gen at2g39030|que hace el gen at2g39030|dime la funcion de at2g39030|explicame at2g39030|cual es el rol de at2g39030)[!., ]*$",
    response:
      "🔬 <strong>AT2G39030</strong> esta moderadamente sobreexpresado (logFC 2.55), y participa en la senalizacion celular y el metabolismo. Su regulacion puede generar efectos en cascada, afectando la expresion de otros genes y modulando rutas de sintesis de proteinas y respuesta a hormonas. Esto sugiere que las celulas pueden estar activando mecanimos de adaptacion o crecimiento en respuesta a la condicion estudiada. Un FDR de 9.39e-48 respalda que esta diferencia es estadisticamente muy confiable.",
  },
  {
    pattern:
      "^(que significa un gen diferencialmente expresado|explicame un gen diferencialmente expresado|para que sirve un gen diferencialmente expresado|que es un gen deg)[!., ]*$",
    response:
      "📊 Un gen diferencialmente expresado (DEG) es aquel que muestra cambios significativos de expresion entre dos condiciones o tratamientos. Detectar DEGs permite identificar genes que responden a perturbaciones, enfermedades o condiciones experimentales, revelando rutas biologicas afectadas y posibles mecanismos celulares alterados. Analizar DEGs ayuda a interpretar como se adaptan las celulas y tejidos, y proporciona una base para estudios funcionales posteriores o potenciales estrategias de intervencion.",
  },
  {
    pattern:
      "^(explicame logfc|que es logfc|para que sirve logfc|como interpretar logfc)[!., ]*$",
    response:
      "📈 El <strong>logFC</strong> (log fold change) indica la magnitud y direccion del cambio en expresion de un gen entre dos condiciones. Valores positivos significan sobreexpresion, valores negativos subexpresion. Por ejemplo, un logFC de -4.16 para AT5G23820 indica una disminucion de 16 veces, mientras que un logFC de 4.2 para AT5G59310 indica un aumento de 18 veces. LogFC combinado con FDR permite evaluar la relevancia biologica y estadistica de cada gen, y se visualiza facilmente en volcano plots y heatmaps.",
  },
  {
    pattern:
      "^(que significa un fdr pequeno|explicame fdr|como interpretar fdr|para que sirve fdr)[!., ]*$",
    response:
      "✅ El <strong>FDR</strong> (false discovery rate) mide la probabilidad de que un gen identificado como diferencialmente expresado sea un falso positivo. Un FDR muy bajo, como 2.9e-123 para AT5G23820, indica que el cambio es altamente confiable. Junto con logFC, permite seleccionar DEGs relevantes y priorizar genes para estudios posteriores, asegurando que los resultados sean reproducibles y biologicamente significativos.",
  },
  {
    pattern:
      "^(que es un volcano plot|explicame volcano plot|para que sirve un volcano plot)[!., ]*$",
    response:
      "🌋 Un <strong>volcano plot</strong> combina la magnitud del cambio (logFC) y la significancia estadistica (FDR o p-value) de cada gen. Genes con cambios grandes y significativos aparecen en los extremos laterales superiores, destacando DEGs relevantes. Por ejemplo, AT5G23820 estaria en el extremo izquierdo y AT5G59310 en el derecho. Este tipo de grafico permite identificar rapidamente genes candidatos y visualizar patrones generales de regulacion en el experimento.",
  },
  {
    pattern:
      "^(que es un heatmap|explicame heatmap|como interpretar heatmap|)[!., ]*$",
    response:
      "🗺️ Un <strong>heatmap</strong> representa la expresion de varios genes en diferentes muestras. Las filas son genes, las columnas muestras, y los colores representan niveles de expresion, desde bajos (colores frios) hasta altos (colores calidos). Los patrones visibles permiten detectar co-regulacion de genes, clusters de muestras similares y tendencias globales de expresion. Por ejemplo, AT5G23820 se veria en colores frios por su baja expresion, mientras que AT5G59310 se veria en colores calidos.",
  },
  {
    pattern:
      "^(como afecta la regulacion de un gen a la celula|explicame efecto de regulacion de un gen|para que sirve regular un gen|que pasa si se regula un gen)[!., ]*$",
    response:
      "⚡ La regulacion de un gen puede tener efectos directos e indirectos sobre funciones celulares. La sobreexpresion puede activar rutas metabolicas, respuesta al estres o crecimiento celular, mientras que la subexpresion puede ralentizar procesos vitales. Por ejemplo, AT5G23820 disminuido podria reducir division celular y crecimiento, mientras que AT5G59310 sobreexpresado podria potenciar rutas de adaptacion metabolica. Esto permite interpretar la biologia del sistema a nivel molecular y celular de manera integrada.",
  },
  {
    pattern:
      "^(por que estudiar genes diferencialmente expresados|importancia de genes deg|para que sirven los genes deg|para que analizar genes deg)[!., ]*$",
    response:
      "🔍 Analizar genes diferencialmente expresados nos proporciona informacion clave sobre los mecanismos biologicos que cambian entre condiciones. Permite identificar genes y rutas criticas, comprender respuestas celulares y desarrollar hipotesis sobre procesos subyacentes. Esta informacion es util tanto para estudios basicos de biologia como para aplicaciones practicas, como identificacion de biomarcadores o targets de intervencion, mostrando el poder de la integracion de estadistica, bioinformatica y biologia molecular.",
  },
  {
    pattern:
      "^(explicame la importancia de los experimentos de transcriptomica|para que sirven los experimentos de transcriptomica|que es la transcriptomica)[!., ]*$",
    response:
      "🧬 Los experimentos de transcriptomica permiten medir la expresion de todos los genes de una celula o tejido en un momento determinado, proporcionando una imagen completa de la actividad biologica. Esto es fundamental porque los genes no actuan de manera aislada: sus niveles de expresion reflejan la respuesta de la celula a cambios ambientales, tratamientos o enfermedades. Analizar patrones de expresion global, identificar genes diferencialmente expresados y agruparlos en rutas funcionales permite entender procesos complejos como desarrollo, diferenciacion celular o adaptacion metabolica. Ademas, estos datos pueden integrarse con proteomica, metabolomica o datos clinicos para generar modelos predictivos, identificar biomarcadores y guiar estrategias de investigacion.",
  },
  {
    pattern:
      "^(como se interpretan los resultados de un analisis rnaseq|explicame interpretacion de resultados rnaseq|que se aprende de un analisis rnaseq)[!., ]*$",
    response:
      "📊 Interpretar un analisis RNA-seq implica varios niveles de informacion. Primero, se observan las diferencias en la expresion de genes entre condiciones usando medidas como logFC y significancia estadistica (p-values o FDR). Luego, se visualizan patrones globales con heatmaps, volcano plots y analisis de componentes principales, lo que permite identificar clusters de genes y muestras similares. Finalmente, se realiza anotacion funcional y analisis de enriquecimiento, conectando genes con procesos biologicos conocidos. Este enfoque multi-nivel transforma datos crudos en conocimiento biologico accionable, permitiendo entender mecanismos celulares, descubrir genes clave y generar hipotesis sobre respuestas a tratamientos o perturbaciones.",
  },
  {
    pattern:
      "^(que es un volcano plot|que es un volcano plot y para que sirve|explicame volcano plot|)[!., ]*$",
    response:
      "🌋 Un volcano plot es una herramienta de visualizacion que combina la magnitud del cambio en expresion (logFC) con la significancia estadistica (p-value o FDR) de cada gen. En el grafico, los genes con cambios grandes y estadisticamente relevantes aparecen en los extremos superiores izquierdo y derecho, permitiendo identificar rapidamente candidatos biologicos interesantes. Este tipo de visualizacion es especialmente util en estudios de RNA-seq o proteomica para priorizar genes o proteinas para estudios funcionales, detectar patrones globales de regulacion y comunicar resultados complejos de forma clara y visual.",
  },
  {
    pattern:
      "^(como se puede usar la bioinformatica para predecir efectos celulares|explicame bioinformatica para predecir efectos celulares|para que sirve la bioinformatica en biologia)[!., ]*$",
    response:
      "🔬 La bioinformatica permite integrar datos de expresion genica, interacciones proteicas, rutas metabolicas y redes de senalizacion para modelar y predecir efectos celulares. Por ejemplo, si varios genes de una ruta metabolica muestran patrones coordinados de expresion, se puede inferir que la actividad de esa ruta esta alterada, lo que podria afectar crecimiento, division celular o respuesta al estres. Ademas, algoritmos de machine learning y redes neuronales pueden aprender patrones complejos de co-expresion, permitiendo predecir efectos combinados de multiples genes o condiciones. Esto convierte grandes datasets en conocimiento biologico accionable, acelerando descubrimientos y optimizando experimentos futuros.",
  },
  {
    pattern:
      "^(por que es importante integrar estadistica y biologia en el analisis de datos|explicame importancia de estadistica en biologia|para que sirve integrar estadistica y biologia)[!., ]*$",
    response:
      "📈 Integrar estadistica y biologia es fundamental para que los resultados sean robustos y biologicamente interpretables. La estadistica permite cuantificar cambios de expresion, controlar errores y evaluar la significancia de los hallazgos, mientras que la biologia proporciona el contexto necesario para comprender el impacto de esos cambios en rutas, procesos celulares y fenotipos. Sin estadistica, cualquier observacion podria ser un falso positivo; sin biologia, los numeros carecen de sentido. La combinacion de ambas disciplinas permite transformar datos masivos en conocimiento util, generar hipotesis solidas y tomar decisiones fundamentadas en investigaciones de transcriptomica, proteomica y sistemas biologicos complejos.",
  },
  {
    pattern:
      "^(que es un gen sobreexpresado|explicame gen sobreexpresado|para que sirve un gen sobreexpresado)[!., ]*$",
    response:
      "📈 Un gen sobreexpresado es aquel cuya expresion aumenta significativamente en una condicion en comparacion con otra. Esta sobreexpresion puede indicar que la celula esta activando rutas metabolicas, de senalizacion o respuestas especificas a un estimulo. Identificar genes sobreexpresados permite comprender procesos biologicos clave, priorizar genes para estudios funcionales y generar hipotesis sobre adaptaciones celulares o posibles mecanismos de regulacion.",
  },
  {
    pattern:
      "^(que es un gen inhibido|explicame gen inhibido|para que sirve un gen inhibido)[!., ]*$",
    response:
      "📉 Un gen inhibido es aquel cuya expresion disminuye significativamente en una condicion frente a otra. Esta disminucion puede reflejar que la celula esta reprimiendo ciertas rutas metabolicas, senalizacion o funciones especificas. Identificar genes inhibidos permite entender como la celula regula negativamente ciertos procesos, detectar rutas que se encuentran menos activas y guiar investigaciones sobre mecanismos de control y adaptacion.",
  },
  {
    pattern:
      "^(por que es importante estudiar genes sobreexpresados o inhibidos|importancia de genes sobreexpresados e inhibidos|para que sirve estudiar expresion diferencial)[!., ]*$",
    response:
      "🔍 Estudiar genes sobreexpresados o inhibidos es clave para entender la biologia subyacente de un sistema. Permite identificar rutas metabolicas, redes de senalizacion y mecanismos de adaptacion a condiciones o tratamientos. Ademas, proporciona informacion valiosa para generar hipotesis sobre el efecto de perturbaciones, descubrir genes candidatos para intervenciones y encontrar posibles biomarcadores en biomedicina, agricultura o investigacion basica.",
  },
  {
    pattern:
      "^(que se puede hacer una vez detectados genes diferencialmente expresados|como analizar genes diferencialmente expresados|para que sirve detectar DEGs)[!., ]*$",
    response:
      "⚡ Una vez detectados genes diferencialmente expresados, se pueden realizar varias acciones: primero, explorar patrones globales de expresion mediante heatmaps y volcano plots; luego, realizar anotacion funcional y analisis de enriquecimiento de rutas para entender su papel biologico; y finalmente, priorizar genes para estudios funcionales, validacion experimental o diseño de intervenciones. Este flujo transforma datos de RNA-seq en conocimiento biologico accionable.",
  },
  {
    pattern:
      "^(como interpretar un gen sobreexpresado o inhibido en un volcano plot|explicame interpretacion de volcano plot para DEGs|como leer genes sobreexpresados e inhibidos)[!., ]*$",
    response:
      "🌋 En un volcano plot, los genes sobreexpresados aparecen normalmente en el extremo derecho, mientras que los genes inhibidos se encuentran en el extremo izquierdo. La altura de los puntos indica la significancia estadistica (p-value o FDR). Esta visualizacion permite identificar rapidamente genes que muestran cambios biologicamente relevantes y estadisticamente confiables, facilitando la seleccion de DEGs para estudios posteriores y la interpretacion de patrones globales de regulacion genetica.",
  },
  {
    pattern:
      "^(que son los metadatos en un experimento rnaseq|explicame metadatos rnaseq|para que sirven los metadatos)[!., ]*$",
    response:
      "📂 Los metadatos son informacion descriptiva sobre cada muestra del experimento, como su nombre, tipo o condiciones experimentales. Permiten organizar, etiquetar y contextualizar los datos crudos, garantizando que las comparaciones entre grupos biologicos sean correctas. Sin metadatos claros, seria dificil interpretar resultados de expresion diferencial, reproducir analisis o generar visualizaciones precisas como heatmaps y volcano plots.",
  },
  {
    pattern:
      "^(que representa la columna filename|para que sirve filename en los metadatos|explicame filename)[!., ]*$",
    response:
      "🗃️ La columna Filename indica el nombre del archivo original de cada muestra, tal como se almacena en el sistema de secuenciacion o en el pipeline de analisis. Esto permite vincular los datos crudos con la informacion descriptiva y asegura que cada lectura o muestra se asocie correctamente con su grupo biologico y su condicion experimental.",
  },
  {
    pattern:
      "^(que representa la columna name en los metadatos|para que sirve name en un experimento rnaseq|explicame name)[!., ]*$",
    response:
      "🏷️ La columna Name proporciona una etiqueta mas amigable o abreviada para cada muestra, facilitando la referenciacion en matrices de expresion, graficos y tablas de resultados. Por ejemplo, PX_0_1 indica la condicion experimental, el tratamiento o genotipo, y la replica, lo que permite organizar la informacion de manera clara y consistente en el analisis.",
  },
  {
    pattern:
      "^(que representa la columna type en los metadatos|para que sirve type en un experimento rnaseq|explicame type)[!., ]*$",
    response:
      "🔬 La columna Type identifica el grupo biologico o genotipo de cada muestra, como mutante o control. Esto es fundamental para realizar comparaciones diferenciales entre condiciones, agrupar replicados biologicos y entender el contexto biologico de los resultados. Sin esta informacion, seria imposible determinar que genes estan regulados en respuesta a un tratamiento o diferencia genetica.",
  },
  {
    pattern:
      "^(como se usan los metadatos para analisis de expresion diferencial|para que sirven los metadatos en DEGs|explicame uso de metadatos en analisis)[!., ]*$",
    response:
      "⚡ Los metadatos permiten definir las comparaciones correctas entre grupos de muestras para identificar genes diferencialmente expresados. Gracias a ellos, podemos agrupar replicados biologicos, separar controles de mutantes, y etiquetar correctamente los resultados en tablas y graficos. Ademas, facilitan la reproducibilidad del analisis y la interpretacion biologica, asegurando que los patrones observados reflejen cambios reales y no errores de organizacion de datos.",
  },
  {
    pattern:
      "^(que es un gen|explicame que es un gen|para que sirve un gen)[!., ]*$",
    response:
      "🧬 Un gen es una unidad de informacion hereditaria que contiene la secuencia de ADN necesaria para producir una molecula funcional, generalmente una proteina o un ARN funcional. Los genes determinan caracteristicas biologicas y fisiologicas, y su expresion puede variar segun condiciones ambientales o geneticas. Estudiar genes permite entender mecanismos celulares, rutas metabolicas y regulacion genetica.",
  },
  {
    pattern: "^(que es el adn|explicame adn|para que sirve el adn)[!., ]*$",
    response:
      "🧪 El ADN (acido desoxirribonucleico) es la molecula que almacena la informacion genetica de los organismos. Contiene las instrucciones necesarias para construir y mantener una celula, y se organiza en cromosomas. Su secuencia determina la identidad genetica y regula la expresion de genes, siendo la base de la herencia y de la diversificacion biologica.",
  },
  {
    pattern: "^(que es el arn|explicame arn|para que sirve el arn)[!., ]*$",
    response:
      "🧫 El ARN (acido ribonucleico) es una molecula que copia temporalmente la informacion del ADN para sintetizar proteinas o cumplir funciones reguladoras en la celula. Existen distintos tipos de ARN, como el mensajero (ARNm), que lleva la informacion genetica, el ribosomico (ARNr), que forma parte de los ribosomas, y el de transferencia (ARNt), que ayuda en la traduccion de proteinas. El ARN es fundamental para la expresion genetica y para muchos procesos celulares.",
  },
  {
    pattern: "^(que es un avs|explicame avs|para que sirve un avs)[!., ]*$",
    response:
      "🧬 Un AVS (allele variant site) es un sitio en el genoma donde existe una variacion entre individuos, como un SNP o mutacion puntual. Estas variaciones pueden afectar la funcion de genes y proteinas, y son importantes para estudiar diversidad genetica, susceptibilidad a enfermedades y respuesta a tratamientos. Identificar AVS permite correlacionar variantes geneticas con fenotipos y procesos biologicos.",
  },
  {
    pattern: "^(que es un otu|explicame otu|para que sirve un otu)[!., ]*$",
    response:
      "🔬 Un OTU (Operational Taxonomic Unit) es una unidad utilizada en metagenomica y estudios de diversidad microbiana para agrupar secuencias similares de ADN o ARN, representando especies o grupos taxonomicos. Los OTUs permiten analizar la composicion y diversidad de comunidades microbianas sin necesidad de identificar cada especie individualmente, facilitando comparaciones entre muestras y condiciones ambientales.",
  },
  {
    pattern:
      "^(que es una proteina|explicame proteina|para que sirve una proteina)[!., ]*$",
    response:
      "🧪 Una proteina es una macromolecula formada por cadenas de aminoacidos que realiza funciones esenciales en la celula, como catalizar reacciones quimicas (enzimas), transportar moleculas, regular procesos celulares o formar estructuras celulares. La informacion para sintetizar proteinas proviene de los genes y su estudio permite entender la biologia molecular y los mecanismos de muchas enfermedades.",
  },
  {
    pattern:
      "^(como interpretar un volcano plot|explicame interpretacion de volcano plot|para que sirve un volcano plot)[!., ]*$",
    response:
      "🌋 Un volcano plot es una herramienta visual que combina dos dimensiones clave: la magnitud del cambio en expresion (logFC) y la significancia estadística (p-value o FDR). Cada punto representa un gen. Los genes con cambios grandes y significativos aparecen en los extremos superior izquierdo (genes inhibidos) y derecho (genes sobreexpresados). Para obtener conclusiones, primero identificamos los genes que cruzan los umbrales de significancia y cambio de expresión; estos son los candidatos principales para estudios funcionales. Luego, podemos relacionarlos con rutas metabólicas o procesos biológicos para entender qué funciones celulares están afectadas. El volcano plot también permite ver patrones globales de regulación: si muchos genes caen hacia un lado, puede indicar una activación o represión coordinada. En combinación con otras visualizaciones como heatmaps o análisis de enriquecimiento, el volcano plot ayuda a priorizar genes para validación experimental y genera hipótesis sobre mecanismos celulares.",
  },
  {
    pattern:
      "^(como interpretar un heatmap|explicame interpretacion de heatmap|para que sirve un heatmap)[!., ]*$",
    response:
      "🗺️ Un heatmap muestra la expresión de genes en distintas muestras usando colores, donde filas son genes y columnas son muestras. La intensidad del color indica el nivel de expresión. Para extraer conclusiones, primero observamos patrones de co-expresión mediante clustering jerárquico: genes que se regulan de manera similar suelen agruparse, así como muestras con perfiles similares. Esto permite identificar subgrupos de genes y de muestras que responden de manera coordinada a un tratamiento o condición genética. Posteriormente, estos clusters se pueden analizar funcionalmente, relacionando genes con rutas o procesos biológicos. Además, los heatmaps permiten detectar anomalías, replicados atípicos o efectos técnicos. Combinando esta información con los resultados de volcano plots, podemos priorizar genes para experimentos de validación y entender mejor la respuesta biológica global.",
  },
  {
    pattern:
      "^(que significa que (el )?gen at5g23820 este subexpresad(o|a)|que significa que at5g23820 este inhibid(o|a)|explicame si at5g23820 esta subexpresad(o|a)|que implica la subexpresion de at5g23820)[!., ]*$",
    response:
      "📉 Cuando el gen AT5G23820 aparece subexpresado, significa que produce menos ARN en una condición frente a otra. Esto se interpreta así:  \n\n1) Evidencia estadística: un logFC negativo indica reducción de expresión y un FDR bajo confirma que el cambio es fiable. Ambos deben cumplirse para considerar la subexpresión como real.  \n\n2) Significado biológico: menos ARN suele reflejar represión transcripcional, menor estabilidad del ARN o pérdida de activadores. Las consecuencias dependen de la función del gen; por ejemplo, si participa en crecimiento, podría disminuir la actividad de esa vía; si actúa en defensa, podría reducir la respuesta inmune.  \n\n3) Impacto en vías: la bajada de un gen puede afectar procesos conectados. Es útil revisar análisis de enriquecimiento y redes para ver qué rutas podrían estar modulándose.  \n\n4) Comprobaciones técnicas: revisar conteos crudos, consistencia entre réplicas, ausencia de batch effects y la posición del gen en volcano plot y heatmap. Genes con conteos muy bajos requieren cautela.  \n\n5) Validación: la subexpresión debe confirmarse con qPCR y, si procede, con análisis proteico o ensayos funcionales como mutantes o silenciamiento.  \n\n6) Conclusión breve para presentar: AT5G23820 muestra una reducción significativa de expresión en la condición X respecto a Y, consistente entre réplicas y respaldada por logFC negativo y FDR bajo. Esto sugiere una modulación de la vía Z y requiere validación adicional.  \n\n7) Precaución: la subexpresión indica asociación, no causalidad. Se debe integrar con otros análisis antes de atribuir un mecanismo concreto.",
  },
  {
    pattern:
      "^(que significa que (el )?gen at5g59310 este sobreexpresad(o|a)|explicame si at5g59310 esta sobreexpresad(o|a)|que implica la sobreexpresion de at5g59310)[!., ]*$",
    response:
      "📈 Cuando el gen AT5G59310 aparece sobreexpresado en un análisis transcriptómico, significa que sus niveles de ARN son más altos en una condición respecto a otra. Esta observación se interpreta de la siguiente manera:  \n\n1) Evidencia estadística: un logFC positivo indica un aumento de expresión y un FDR bajo confirma que el cambio es robusto y poco probable por azar. La combinación de ambos factores permite afirmar que la sobreexpresión es real.  \n\n2) Significado biológico: más ARN puede reflejar activación transcripcional, mayor estabilidad del ARN, presencia de reguladores activadores o una respuesta adaptativa al estímulo analizado. Las implicaciones dependen de la función del gen. Por ejemplo, si AT5G59310 participa en procesos metabólicos o de defensa, su aumento podría potenciar esas rutas.  \n\n3) Impacto en vías y redes: la sobreexpresión de un gen puede amplificar rutas completas o modificar la actividad de redes de regulación. Es recomendable revisar análisis de enriquecimiento funcional, coexpresión e interacciones para identificar qué procesos podrían estar intensificándose.  \n\n4) Verificación técnica: antes de concluir, conviene revisar conteos crudos, consistencia entre réplicas, calidad de la normalización y la posición del gen en volcano plot y heatmap. Genes con conteos bajos deben interpretarse con prudencia.  \n\n5) Validación experimental: una sobreexpresión detectada computacionalmente debe confirmarse mediante qPCR y, si es relevante, mediante análisis proteico. Ensayos funcionales como líneas transgénicas de sobreexpresión o silenciamiento ayudan a determinar si el aumento tiene consecuencias fenotípicas.  \n\n6) Conclusión breve para presentar: AT5G59310 muestra un incremento significativo de expresión en la condición X frente a Y, respaldado por logFC positivo y un FDR fiable. Esto apunta a una activación de la vía Z en la que participa el gen, pendiente de validación adicional.  \n\n7) Precaución: la sobreexpresión indica asociación con la condición experimental, no causalidad directa. Para proponer mecanismos concretos es necesario integrar más evidencia y realizar experimentos funcionales.",
  },
  {
    pattern:
      "^(que es (un )?ma plot|para que sirve un ma plot|explicame el ma plot)[!., ]*$",
    response:
      "📊 El MA Plot es un gráfico que representa la relación entre la expresión media de los genes (eje X, logCPM) y el cambio relativo en su expresión (eje Y, logFC). Su utilidad principal es identificar genes diferencialmente expresados: aquellos alejados de la línea de y = 0 son los que muestran cambios relevantes. Permite detectar patrones generales de expresión, visualizar sesgos dependientes de la abundancia y localizar genes sobreexpresados o inhibidos de manera rápida.",
  },
  {
    pattern:
      "^(como interpreto los genes que estan arriba o abajo en el ma plot|como saber que genes son importantes en el ma plot|dime que genes mirar en el ma plot)[!., ]*$",
    response:
      "🔎 En el MA Plot, los genes que se encuentran lejos de la línea roja (y = 0) y con logFC alto o bajo son los más relevantes. Para identificarlos, fíjate en:  \n1) Genes con FDR < 0.05 y |logFC| > 1 → considerados diferencialmente expresados.  \n2) Genes con logFC extremos → cambios fuertes, aunque haya que revisar su expresión media.  \n3) Genes con logCPM suficiente → descarta cambios en genes muy poco expresados, que podrían ser ruido.  \n4) Genes de interés biológico → aquellos implicados en vías o funciones relevantes para el experimento.",
  },
  {
    pattern:
      "^(que significa que un gen este sobreexpresado en el ma plot|como saber si un gen esta sobreexpresado en el ma plot|explicame la sobreexpresion de un gen en el ma plot)[!., ]*$",
    response:
      "📈 Cuando un gen aparece sobreexpresado en el MA Plot, significa que su nivel de ARN es mayor en una condición respecto a otra. Un logFC positivo indica aumento de expresión, y un FDR bajo confirma que el cambio es estadísticamente significativo. La interpretación biológica depende de la función del gen: puede reflejar activación de una vía, respuesta adaptativa o regulación transcripcional. Se recomienda revisar la expresión media, consistencia entre réplicas y complementar con volcano plot y heatmap para verificar confiabilidad.",
  },
  {
    pattern:
      "^(que significa que un gen este inhibido en el ma plot|como saber si un gen esta inhibido en el ma plot|explicame la inhibicion de un gen en el ma plot)[!., ]*$",
    response:
      "📉 Cuando un gen aparece inhibido en el MA Plot, significa que su nivel de expresión es menor en una condición frente a otra. Un logFC negativo indica disminución de expresión, y un FDR bajo asegura que el cambio es estadísticamente significativo. Biológicamente, puede indicar represión transcripcional, inhibición de vías o regulación por factores represores. Como siempre, conviene revisar la expresión media, replicados y visualizaciones complementarias para confirmar la observación.",
  },
  {
    pattern:
      "^(como identificar genes con cambios significativos en el ma plot|dime que genes son relevantes en el ma plot|que criterios usar para seleccionar genes del ma plot)[!., ]*$",
    response:
      "✅ Para identificar genes relevantes en el MA Plot, aplica los siguientes criterios:  \n1) FDR < 0.05 → garantiza significancia estadística.  \n2) |logFC| > 1 → asegura un cambio de expresión considerable.  \n3) logCPM suficiente → descarta genes con muy baja expresión.  \n4) Prioridad a genes de interés biológico → relacionados con vías o funciones esperadas en el experimento.  \nLos genes que cumplen estos criterios suelen aparecer lejos de la línea y = 0 y son los que conviene explorar primero en análisis funcionales o experimentales.",
  },

  // PREGUNTAS EDU
  {
    pattern:
      "^(que es un fichero fastq|que es un archivo fastq|explicame que es fastq|que significa fastq|para que sirve un fastq|definicion de fastq)[!., ]*$",
    response:
      "📁 Un fichero FASTQ es un formato fundamental en secuenciacion masiva que almacena dos piezas clave de informacion: la secuencia de nucleotidos obtenida por el secuenciador y una linea paralela con las calidades asociadas a cada base. Cada lectura se representa en cuatro lineas: identificador, secuencia, separador y cadena de calidades codificadas en ASCII. Este formato permite evaluar no solo el contenido biologico de las lecturas sino tambien la confianza tecnica en cada posicion. Los archivos FASTQ son el punto de partida de cualquier pipeline de NGS y sirven para control de calidad, alineamiento, cuantificacion y posteriores analisis bioinformaticos.",
  },
  {
    pattern:
      "^(que es una lectura en ngs|que es una lectura|que significa una lectura|explicame que es una lectura en secuenciacion|definicion de lectura ngs)[!., ]*$",
    response:
      "🧬 Una lectura en NGS es una pequena secuencia de nucleotidos generada por la maquina de secuenciacion a partir de un fragmento de adn o arn. Cada lectura representa una observacion directa del material genetico fragmentado y es la unidad basica de informacion con la que trabaja un pipeline de analisis. Las lecturas pueden ser de un extremo (single-end) o de dos extremos (paired-end), lo que permite reconstrucciones mas precisas. Su calidad, longitud y cantidad determinan la fiabilidad del alineamiento, la cuantificacion y la interpretacion biologica final.",
  },
  {
    pattern:
      "^(que es un adaptador en ngs|que es un adaptador|que significa adaptador en secuenciacion|explicame que es un adaptador de libreria|para que sirven los adaptadores ngs)[!., ]*$",
    response:
      "🧷 Un adaptador es una secuencia sintetica que se anade a los extremos de los fragmentos de adn o arn durante la preparacion de librerias para secuenciacion masiva. Estos adaptadores permiten que los fragmentos se unan al flujo de celulas del secuenciador, posibilitan su amplificacion y contienen elementos necesarios para el reconocimiento por parte de la plataforma. Aunque esenciales, a veces aparecen dentro de las lecturas si el fragmento biologico es muy corto; por eso se realiza un recorte de adaptadores para evitar errores en los pasos posteriores como el alineamiento o la cuantificacion.",
  },
  {
    pattern:
      "^(que es rnaseq|que es rna seq|que es rna-seq|explicame que es rnaseq|definicion de rnaseq|para que sirve rnaseq)[!., ]*$",
    response:
      "🔬 RNA-seq es una tecnica de secuenciacion masiva que permite detectar, cuantificar y caracterizar el arn presente en una muestra biologica. Mediante la conversion del arn en cdna y su posterior secuenciacion, es posible medir niveles de expresion genica, identificar isoformas, descubrir transcritos nuevos y analizar cambios en la actividad genetica entre condiciones. RNA-seq es una herramienta central en transcriptomica porque ofrece una vision global, cuantitativa y altamente sensible del paisaje de expresion celular.",
  },
  {
    pattern:
      "^(que es una proteina|explicame que es una proteina|definicion de proteina|para que sirve una proteina|que significa proteina)[!., ]*$",
    response:
      "⚙️ Una proteina es una molecula funcional formada por cadenas de aminoacidos que adopta estructuras tridimensionales especificas para desempenar tareas esenciales en la celula. Las proteinas actuan como enzimas, receptores, factores estructurales, reguladores, transportadores y motores moleculares, entre muchos otros roles. Su secuencia esta codificada por los genes, de modo que los cambios en la expresion genica, como los estudiados mediante rna-seq, suelen traducirse en variaciones en la actividad proteica y en los procesos celulares que estas regulan.",
  },
  {
    pattern:
      "^(que mide rnaseq|que mide rna seq|que mide rna-seq|que informa rnaseq|para que sirve rnaseq medir|que cuantifica rnaseq)[!., ]*$",
    response:
      "📊 RNA-Seq mide la expresión génica cuantificando los niveles de ARN transcrito en una muestra. Para ello, convierte las moléculas de ARN en cDNA, las fragmenta y las secuencia, generando lecturas que representan la actividad de cada gen. A partir de estas lecturas es posible comparar condiciones, detectar genes diferencialmente expresados, identificar isoformas y evaluar la actividad global de las rutas biológicas.",
  },
  {
    pattern:
      "^(que es un gen|definicion de gen|explicame que es un gen|que significa gen|como se define un gen)[!., ]*$",
    response:
      "🧬 Un gen es una región del ADN que contiene la información necesaria para producir una proteína o una molécula de ARN funcional. Los genes incluyen regiones codificantes y elementos reguladores que controlan cuándo y cuánto se expresa. Su actividad determina procesos celulares esenciales como crecimiento, señalización, metabolismo y respuesta a estrés.",
  },
  {
    pattern:
      "^(que es un transcrito|que significa transcrito|definicion de transcrito|que es un arn transcrito|explicame que es un transcrito)[!., ]*$",
    response:
      "📄 Un transcrito es una molécula de ARN generada a partir de un gen mediante el proceso de transcripción. Puede corresponder a ARN mensajero, ARN no codificante o isoformas alternativas de un mismo gen. Los transcritos reflejan la actividad real del gen en un momento y condición determinados, por lo que su cuantificación es fundamental en estudios de expresión.",
  },
  {
    pattern:
      "^(que es una libreria de rnaseq|que es la libreria rnaseq|explicame libreria de rnaseq|como se forma una libreria rnaseq|definicion de libreria de secuenciacion)[!., ]*$",
    response:
      "📚 Una librería de RNA-Seq es el conjunto de fragmentos de ARN convertidos en cDNA, adaptados y preparados para ser secuenciados. Durante su construcción se fragmenta el ARN, se añaden adaptadores, se seleccionan tamaños y se amplifica el material. La calidad de la librería determina directamente la profundidad, diversidad y fiabilidad de todo el experimento de secuenciación.",
  },
  {
    pattern:
      "^(que es (una )?lectura|(que significa|explicame|definicion) (una )?lectura( en ngs)?|que representa (una )?lectura)[!., ]*$",
    response:
      "🧩 Una lectura es una pequeña secuencia de nucleótidos generada por el secuenciador a partir de un fragmento de ADN o ARN. Representa la unidad básica de información de cualquier análisis NGS. Su calidad, longitud y distribución influyen directamente en el alineamiento, la cuantificación y la interpretación biológica del experimento.",
  },
  {
    pattern:
      "^(que es (un )?adaptador|(que significa|explicame|definicion) (un )?adaptador( en ngs)?|para que sirven los adaptadores)[!., ]*$",
    response:
      "🧷 Un adaptador es una secuencia sintética que se añade a los extremos de los fragmentos durante la preparación de librerías. Permite que los fragmentos se unan a la celda de flujo, se amplifiquen y sean reconocidos por la plataforma de secuenciación. Si aparecen dentro de las lecturas, deben recortarse para evitar errores en el alineamiento o en la cuantificación.",
  },
  {
    pattern:
      "^(que es (la )?alineacion|(que significa|explicame) (el )?alineamiento|definicion de alineacion( ngs)?|para que sirve alinear lecturas)[!., ]*$",
    response:
      "🧭 La alineación es el proceso de mapear cada lectura contra un genoma o transcriptoma de referencia para determinar su origen biológico. Este paso permite saber a qué gen o región pertenece cada fragmento secuenciado. Una alineación precisa es crucial para cuantificar correctamente la expresión, detectar variantes y explorar la estructura del transcriptoma.",
  },
  {
    pattern:
      "^(que es (un )?genoma de referencia|(para que sirve|explicame|definicion) (un )?genoma de referencia)[!., ]*$",
    response:
      "🗺️ Un genoma de referencia es una secuencia completa y anotada que actúa como guía para interpretar las lecturas de un experimento de secuenciación. Permite ubicar cada lectura en un contexto biológico, identificar genes, regiones reguladoras y variantes. Su calidad y actualización afectan directamente a la precisión del alineamiento y la anotación.",
  },
  {
    pattern:
      "^(que es (la )?cuantificacion|(que significa|explicame|definicion) (la )?cuantificacion( genica| rnaseq)?|para que sirve (la )?cuantificacion)[!., ]*$",
    response:
      "🔢 La cuantificación consiste en contar cuántas lecturas se asignan a cada gen o transcrito para estimar su nivel de expresión. Este paso convierte datos de secuenciación en valores numéricos comparables entre muestras. La cuantificación puede ser a nivel de genes o de isoformas y suele corregirse mediante normalizaciones para evitar sesgos técnicos.",
  },
  {
    pattern:
      "^(que es (tpm)|que significa (tpm)|definicion de (tpm)|explicame que es (tpm)|para que sirve (tpm) en rnaseq)[!., ]*$",
    response:
      "📐 TPM significa Transcripts Per Million y es una medida de expresión normalizada que ajusta los conteos de lecturas según la longitud del transcrito y la profundidad total de secuenciación. Esto permite comparar niveles de expresión entre genes dentro de una misma muestra y reduce sesgos derivados de diferencias técnicas entre librerías.",
  },
  {
    pattern:
      "^(que es (fpkm)|que significa (fpkm)|que quiere decir (fpkm)|definicion de (fpkm))[!., ]*$",
    response:
      "📐 FPKM (Fragments Per Kilobase per Million) es una medida normalizada de expresión que ajusta el número de fragmentos detectados por la longitud del transcrito (en kilobases) y por la profundidad total de secuenciación (por cada millón de fragmentos). Su objetivo es permitir comparaciones relativas de expresión entre genes en una misma muestra, corrigiendo sesgos debidos a la longitud del gen y a diferentes profundidades de secuenciación. Para comparaciones entre muestras se suelen preferir otras normalizaciones más modernas, como TPM o métodos basados en conteos normalizados.",
  },
  {
    pattern:
      "^(que es (la )?normalizacion|que significa normalizacion en rnaseq|para que sirve (la )?normalizacion|explicame normalizacion de datos)[!., ]*$",
    response:
      "⚖️ La normalización es el proceso de ajustar los datos de secuenciación para eliminar sesgos técnicos que no están relacionados con la biología real, como diferencias en profundidad de secuenciación, tamaños de librería o composición de la muestra. El objetivo es que las diferencias observadas entre muestras reflejen cambios biológicos y no artefactos. Los métodos de normalización incluyen ajustes por tamaño de librería, TMM, rlog, VST o TPM, y la elección depende del análisis posterior y de la naturaleza de los datos.",
  },
  {
    pattern:
      "^(que es (un )?gen diferencialmente expresado|que es (un )?deg|que significa gen diferencialmente expresado|definicion de deg)[!., ]*$",
    response:
      "🔍 Un gen diferencialmente expresado (DEG) es aquel que muestra un cambio estadísticamente significativo en su nivel de expresión entre dos o más condiciones. Para considerarlo DEG se combinan la magnitud del cambio (por ejemplo logFC) y la significancia estadística (p-value y FDR). Identificar DEGs permite señalar genes candidatos implicados en la respuesta biológica al tratamiento, condición o fenotipo estudiado.",
  },
  {
    pattern:
      "^(que es (un )?control biologico|que significa control biologico|para que sirve (un )?control biologico|definicion de control biologico)[!., ]*$",
    response:
      "🔬 Un control biológico es una muestra que actúa como referencia para interpretar efectos experimentales: puede ser el estado no tratado, el tipo salvaje o una condición baseline. Su función es proporcionar un punto de comparación para identificar cambios reales debidos al tratamiento o condición. Controles biológicos bien elegidos aumentan la validez biológica y facilitan la interpretación de DEGs.",
  },
  {
    pattern:
      "^(que es (un )?replicado|que significa (un )?replicado biologico|para que sirven (los )?replicados|definicion de (un )?replicado|explicame replicados en rnaseq)[!., ]*$",
    response:
      "🔁 Un replicado es una repetición biológica o técnica de una medición experimental. Los replicados biológicos capturan la variabilidad natural entre individuos o muestras y son esenciales para evaluar la reproducibilidad y la significancia estadística. Los replicados técnicos evalúan la variabilidad del proceso de medición. En RNA-Seq, disponer de varios replicados biológicos mejora la potencia estadística para detectar DEGs.",
  },
  {
    pattern:
      "^(que es (un )?pca|que significa (pca)|para que sirve (el )?pca|explicame (el )?pca en analisis de datos|que es el analisis de componentes principales)[!., ]*$",
    response:
      "📉 PCA (Análisis de Componentes Principales) es una técnica de reducción de dimensiones que resume la variabilidad de un dataset en unas pocas componentes lineales. En RNA-Seq se usa para explorar la estructura global de los datos, detectar agrupamientos de muestras, identificar outliers y evaluar efectos como tratamiento o batch. Mostrar un PCA claro con replicados agrupados respalda la calidad del experimento y permite interpretar patrones globales de expresión génica.",
  },
  {
    pattern:
      "^(que es (un )?fastqc|que hace (fastqc)|para que sirve (fastqc)|explicame (fastqc)|control de calidad fastq)[!., ]*$",
    response:
      "🛠 FastQC es una herramienta de control de calidad para archivos FASTQ que ofrece reportes sobre la calidad por base, longitud de lecturas, contenido de GC, presencia de adaptadores y otros metadatos. Ayuda a detectar problemas técnicos antes del alineamiento, permitiendo decidir si es necesario realizar trimming, filtrado o repetir la obtención de librerías.",
  },
  {
    pattern:
      "^(que es (el )?trimming|que significa trimming en ngs|para que sirve (el )?trimming|explicame trimming de lecturas|recorte de adaptadores y bases)[!., ]*$",
    response:
      "✂️ Trimming es el proceso de recortar o eliminar porciones no deseadas de las lecturas, típicamente adaptadores y bases de baja calidad en los extremos. Este paso mejora la calidad del alineamiento y reduce falsos emparejamientos. Las herramientas de trimming permiten definir umbrales de calidad y tamaños mínimos para mantener lecturas fiables y listas para el análisis.",
  },
  {
    pattern:
      "^(que es (el )?star|que hace (star)|para que sirve (star)|explicame star alineador|alineador star rnaseq)[!., ]*$",
    response:
      "⚡ STAR es un alineador rápido y sensible para RNA-Seq diseñado para mapear lecturas contra un genoma de referencia, manejando eficientemente empalmes (splicing) y lecturas emparejadas. Se distingue por su velocidad y capacidad para detectar empalmes noveles, siendo una opción ampliamente usada en pipelines de RNA-Seq cuando se requiere mapeo al genoma con alta precisión.",
  },
  {
    pattern:
      "^(que es (el )?hisat2|que hace (el )?hisat2|para que sirve (el )?hisat2|explicame (el )?hisat2 alineador|alineador hisat2 rnaseq)[!., ]*$",
    response:
      "🧭 HISAT2 es un alineador eficiente para RNA-Seq que utiliza índices basados en grafos y estructuras compactas para mapear lecturas con bajo uso de memoria. Es rápido y maneja bien empalmes y variaciones polimórficas, siendo útil cuando se busca equilibrio entre velocidad, memoria y sensibilidad en el mapeo al genoma o transcriptoma.",
  },
  {
    pattern:
      "^(que es (el )?featurecounts|que hace (el )?featurecounts|para que sirve (el )?featurecounts|explicame (el )?featurecounts|cuantificacion de genes featurecounts)[!., ]*$",
    response:
      "🧮 FeatureCounts es una herramienta para cuantificar lecturas asignadas a características génicas (genes, exones) a partir de archivos de alineamiento. Realiza un conteo rápido y eficiente, produciendo matrices de conteos que son la base para análisis diferencial con paquetes como DESeq2 o edgeR.",
  },
  {
    pattern:
      "^(que es (el )?salmon|que hace (el )?salmon|para que sirve (el )?salmon|explicame (el )?salmon cuantificacion|cuantificacion de transcritos salmon)[!., ]*$",
    response:
      "🐟 Salmon es un software para cuantificación rápida y precisa de transcritos que trabaja en modo alineamiento rápido o alignment-free mediante modelos de k-mers y corrección de sesgos. Devuelve estimaciones de abundancia a nivel de transcrito y es muy utilizado por su velocidad y por producir valores normalizados listos para análisis posteriores.",
  },
  {
    pattern:
      "^(que es (el )?deseq2|que hace (el )?deseq2|para que sirve (el )?deseq2|explicame (el )?deseq2|analisis diferencial deseq2)[!., ]*$",
    response:
      "📈 DESeq2 es un paquete en R para el análisis diferencial de expresión basado en modelos de conteo (modelo negativo binomial). Normaliza los conteos, estima dispersiones y calcula logFC y valores p ajustados (FDR), ofreciendo herramientas para visualización y filtrado. Es muy usado por su robustez y fácil integración en pipelines de RNA-Seq.",
  },
  {
    pattern:
      "^(que es (el )?edger|que hace (el )?edger|para que sirve (el )?edger|explicame (el )?edger|analisis diferencial edger)[!., ]*$",
    response:
      "🔎 edgeR es una suite estadística para el análisis diferencial de datos de conteo que también utiliza modelos basados en la distribución binomial negativa. Incluye métodos de normalización, estimación de dispersión y pruebas exactas o por regresión para detectar DEGs, siendo una alternativa complementaria a DESeq2.",
  },
  {
    pattern:
      "^(que es (un )?heatmap|que significa (un )?heatmap|para que sirve (un )?heatmap|explicame (un )?heatmap|interpretacion de heatmap)[!., ]*$",
    response:
      "🗺️ Un heatmap es una visualización que muestra niveles de expresión de muchos genes a través de varias muestras usando una escala de colores. Las filas representan genes, las columnas representan muestras, y el clustering jerárquico de filas y columnas permite identificar patrones de coexpresión y subgrupos de muestras. Los heatmaps son útiles para explorar co-regulación, detectar outliers y seleccionar clusters para análisis funcional.",
  },
  {
    pattern:
      "^(que es (un )?volcano plot|que es (un )?volcanoplot|que es (un )?volcano|explicame (un )?volcano plot|para que sirve (un )?volcano plot|como interpretar (un )?volcano plot)[!., ]*$",
    response:
      "🌋 Un volcano plot es un gráfico que combina la magnitud del cambio en expresión (logFC) y la significancia estadística (p-value o FDR) para cada gen. En el eje X aparece el logFC y en el eje Y la significancia, de modo que los genes con cambios grandes y robustos quedan en los extremos superiores. Para sacar conclusiones se identifican primero los genes que superan umbrales de logFC y FDR, luego se priorizan según relevancia biológica y se integran con análisis de rutas y heatmaps para proponer hipótesis.",
  },
  {
    pattern:
      "^(que es (la )?anotacion funcional|que significa (la )?anotacion funcional|para que sirve (la )?anotacion funcional|explicame (la )?anotacion funcional|definicion de (la )?anotacion funcional)[!., ]*$",
    response:
      "🧩 La anotación funcional es el proceso de asignar términos biológicos, rutas o funciones a genes o transcritos basándose en bases de datos y evidencias. Incluye clasificaciones como ontologías, rutas metabólicas y categorías moleculares, y permite interpretar listas de DEGs en términos biológicos. Con anotación funcional se pasa de una lista de genes a procesos y mecanismos potenciales que explican los cambios observados.",
  },
  {
    pattern:
      "^(que es (el )?go|que significa (el )?go|que es (el )?gene ontology|explicame (el )?go|para que sirve (el )?go|definicion de (el )?gene ontology)[!., ]*$",
    response:
      "📚 GO (Gene Ontology) es un sistema estandarizado que clasifica funciones biológicas, procesos y componentes celulares asociados a genes. GO facilita la anotación y el análisis de enriquecimiento, permitiendo agrupar genes según funciones compartidas y detectar procesos sobre o subrepresentados en un experimento. Es una herramienta clave para interpretar resultados de transcriptómica a nivel funcional.",
  },
  {
    pattern:
      "^(que es (el )?kegg|que significa (el )?kegg|que es (la )?base de datos kegg|explicame (el )?kegg|para que sirve (el )?kegg|definicion de (el )?kegg)[!., ]*$",
    response:
      "🗺 KEGG es una base de datos de rutas metabólicas y de señalización que organiza genes y proteínas en mapas funcionales. Se utiliza para identificar rutas afectadas en un experimento y para visualizar cómo cambios en la expresión pueden alterar procesos metabólicos o de señalización. KEGG complementa la anotación funcional y ayuda a generar hipótesis mecanísticas sobre los efectos observados.",
  },
  {
    pattern:
      "^(que es (el )?batch effect|que significa batch effect|que es (el )?efecto (de )?lote|que es (el )?efecto batch|explicame batch effect|como afecta batch effect|para que sirve identificar batch effect)[!., ]*$",
    response:
      "⚠️ El batch effect es una variación no biológica introducida por diferencias técnicas entre lotes de muestras (fecha, operador, reactivos, secuenciador). Puede enmascarar o confundir señales biológicas si no se identifica y corrige. Para mitigarlo se utilizan diseños experimentales balanceados, controles y métodos estadísticos de corrección como COMBAT o inclusión de covariables en el modelo.",
  },
  {
    pattern:
      "^(que es (el )?filtrado de genes|que significa filtrado de genes|para que sirve filtrar genes|explicame filtrado de genes|definicion de filtrado de genes)[!., ]*$",
    response:
      "🧹 El filtrado de genes consiste en eliminar genes con expresión muy baja o poco informativa antes del análisis diferencial. Esto reduce el ruido, mejora la potencia estadística y evita resultados espurios causados por conteos escasos. Las reglas de filtrado pueden basarse en umbrales de conteos por millón, número mínimo de muestras con expresión o en la distribución global de los datos.",
  },
  {
    pattern:
      "^(que es (el )?rpkm|que significa rpkm|definicion de rpkm|que es reads per kilobase per million|explicame rpkm)[!., ]*$",
    response:
      "📏 RPKM (Reads Per Kilobase per Million) es una medida de expresión que normaliza el número de reads por la longitud del gen y por la profundidad de secuenciación. Permite comparar la expresión relativa de genes dentro de una muestra, pero tiene limitaciones para comparaciones entre muestras frente a métodos más modernos como TPM o normalizaciones basadas en conteos.",
  },
  {
    pattern:
      "^(que es (un )?transcriptoma|que significa transcriptoma|definicion de transcriptoma|explicame transcriptoma|para que sirve transcriptoma)[!., ]*$",
    response:
      "🌐 El transcriptoma es el conjunto completo de transcritos presentes en una célula o tejido en un momento concreto, incluyendo ARN mensajero y ARN no codificante. Estudiar el transcriptoma mediante RNA-Seq ofrece una visión global de la actividad genética, permite identificar isoformas, transcritos nuevos y cambios dinámicos en respuesta a condiciones o tratamientos.",
  },
  {
    pattern:
      "^(que es (un )?polirribosoma|que es polirribosoma o polirradio|que significa polirribosoma|explicame polirribosoma|para que sirve polirribosoma)[!., ]*$",
    response:
      "⚙️ Un polirribosoma es un complejo formado por varios ribosomas unidos al mismo ARN mensajero durante la traducción. La presencia de polirribosomas indica traducción activa y permite inferir el nivel de producción proteica de un transcrito. Técnicas como perfilado de ribosomas (Ribo-Seq) explotan este concepto para estudiar la eficiencia de traducción a nivel génico.",
  },
  {
    pattern:
      "^(que es (la )?ribodeplecion|que significa ribodeplecion|para que sirve (la )?ribodeplecion|explicame ribodeplecion en rnaseq)[!., ]*$",
    response:
      "🧪 La ribodepleción es un paso de preparación de librerías que elimina ARN ribosomal (rRNA) antes de secuenciar, reduciendo la proporción de secuencias ribosomales y aumentando la cobertura de transcritos de interés. Es especialmente útil en muestras con alto contenido de rRNA o cuando se busca estudiar ARN no codificante y transcritos poco abundantes.",
  },
  {
    pattern:
      "^(que es (la )?polya selection|que significa polya selection|explicame polya selection|para que sirve polya selection|que es seleccion polya)[!., ]*$",
    response:
      "🎯 PolyA selection es un método de enriquecimiento de ARN mensajero que captura transcritos con cola poli-A usando oligos de timina unidos a perlas o columnas. Ventajas: aumenta la fracción de ARN codificante en la librería y mejora la eficiencia para estudiar mRNA. Limitaciones: no captura bien ARN sin cola poli-A (p.ej., algunos ARN no codificantes) y puede introducir sesgos hacia los extremos 3' si la integridad del ARN es baja. Alternativa: ribodepleción, útil cuando se quiere estudiar ARN no codificante o muestras degradadas.",
  },
  {
    pattern:
      "^(que es (el )?single end|que es single end sequencing|que significa single end|explicame single end|para que sirve single end|que es secuenciacion single end)[!., ]*$",
    response:
      "🔬 Single-end es la secuenciación en la que se lee un solo extremo de cada fragmento. Ventajas: menor coste y flujo de trabajo más simple. Inconvenientes: menor información sobre la posición precisa y la estructura del fragmento comparado con paired-end, lo que puede dificultar la resolución de isoformas y el mapeo en regiones repetitivas. Se usa cuando el objetivo es cuantificación básica y ahorro de costes.",
  },
  {
    pattern:
      "^(que es (el )?paired end|que es paired end sequencing|que significa paired end|explicame paired end|para que sirve paired end|que es secuenciacion paired end)[!., ]*$",
    response:
      "🔗 Paired-end consiste en secuenciar ambos extremos de cada fragmento, obteniendo dos lecturas por fragmento. Esto mejora el mapeo en regiones repetitivas, facilita la detección de empalmes y isoformas, y permite estimar el tamaño del fragmento (insert size). Es la opción preferida cuando se busca mayor resolución estructural o para ensamblados y detección de variantes.",
  },
  {
    pattern:
      "^(que es (la )?longitud de lectura|que significa longitud de lectura|que es read length|explicame longitud de lectura|para que sirve longitud de lectura|cuantas bases tiene una lectura)[!., ]*$",
    response:
      "📏 La longitud de lectura es el número de nucleótidos en cada lectura generada por el secuenciador. Lecturas más largas mejoran el mapeo unívoco y la resolución de isoformas, pero pueden tener mayor tasa de error en las posiciones finales. La elección del read length depende del objetivo: estudios de cuantificación básica toleran lecturas cortas, mientras que análisis de isoformas o ensamblaje se benefician de lecturas largas.",
  },
  {
    pattern:
      "^(que es (la )?profundidad de secuenciacion|que significa profundidad de secuenciacion|que es sequencing depth|explicame profundidad de secuenciacion|para que sirve profundidad de secuenciacion|que significa cobertura de secuenciacion)[!., ]*$",
    response:
      "🔋 La profundidad de secuenciación (o cobertura) es la cantidad total de lecturas generadas por muestra. Mayor profundidad aumenta la sensibilidad para detectar transcritos raros y mejora la cuantificación, pero conlleva mayores costes. La profundidad necesaria depende del objetivo: estudios de genes muy expresados requieren menos reads; detección de transcritos poco abundantes, isoformas o empalmes requiere mayor profundidad. Conviene evaluar curvas de saturación para decidir el punto eficiente.",
  },
  {
    pattern:
      "^(que es (el )?multi mapping|que significa multi mapping|que es multimapping|explicame multi mapping|como afectan las lecturas multimapeadas|que son lecturas multimapeadas)[!., ]*$",
    response:
      "🌀 Multi-mapping ocurre cuando una lectura puede alinearse en varias posiciones del genoma o transcriptoma por similitud de secuencia (regiones repetidas, familias de genes, isoformas). Esto complica la asignación de lecturas a genes y la cuantificación. Estrategias: descartar lecturas multimapeadas, asignarlas probabilísticamente (herramientas como Salmon) o usar mapeo con parámetros que minimicen multi-mapping. Es importante reportar y manejar estas lecturas para evitar sesgos.",
  },
  {
    pattern:
      "^(que es (un )?splice junction|que es splice junction|que significa splice junction|explicame splice junction|que es una union de empalme|como se detecta splice junction)[!., ]*$",
    response:
      "✂️ Un splice junction es la posición donde dos exones quedan unidos en el ARN maduro tras el empalme (splicing). La detección de splice junctions en RNA-Seq permite identificar isoformas, empalmes alternativos y empalmes noveles no anotados. Herramientas de alineamiento sensibilizadas en splicing (p.ej., STAR, HISAT2) reportan lecturas que atraviesan junctions, clave para el estudio de regulación alternativa del splicing.",
  },
  {
    pattern:
      "^(que es (un )?intron|que significa intron|que es intron en arn|explicame intron|definicion de intron|como se define intron)[!., ]*$",
    response:
      "📦 Un intrón es una región no codificante dentro de un gen que se transcribe en pre-mRNA pero se elimina durante el procesamiento para generar el ARN maduro. La presencia de lecturas en intrones puede indicar transcritos nucleares, retención de intrones o contaminación por ARN no procesado. Interpretar lecturas intrónicas requiere considerar el protocolo de preparación (PolyA selection vs ribodepleción) y la integridad del ARN.",
  },
  {
    pattern:
      "^(que es (un )?exon|que significa exon|que es exon en arn|explicame exon|definicion de exon|como se define exon)[!., ]*$",
    response:
      "🔹 Un exón es una región del gen que permanece en el ARN maduro y que suele codificar aminoácidos en el caso de mRNA. Las lecturas que mapean a exones se usan para cuantificar la expresión génica, y las lecturas que atraviesan junctions entre exones ayudan a definir isoformas. La anotación exónica es esencial para asignar conteos correctamente y para el análisis de isoformas.",
  },
  {
    pattern:
      "^(que es (una )?isoforma|que significa isoforma|que es isoforma de un gen|explicame isoforma|definicion de isoforma|como se define isoforma)[!., ]*$",
    response:
      "🔁 Una isoforma es una variante de transcrito originada por eventos como empalme alternativo, uso alternativo de promotores o terminadores, que produce transcritos con exones distintos o en distinto orden. Las isoformas pueden tener diferencias funcionales, subcelulares o de estabilidad. Cuantificar isoformas es más complejo que cuantificar genes y requiere herramientas y lecturas (preferiblemente paired-end y con suficiente longitud) que permitan discriminar variantes.",
  },
  {
    pattern:
      "^(que es (la )?cuantificacion pseudoalineada|que es cuantificacion pseudoalineada|explicame cuantificacion pseudoalineada|para que sirve (la )?cuantificacion pseudoalineada|que significa pseudoalineacion)[!., ]*$",
    response:
      "🚀 La cuantificación pseudoalineada es un enfoque rápido para estimar la abundancia de transcritos sin realizar un alineamiento completo al genoma. Herramientas como Salmon o Kallisto comparan k-mers de las lecturas con el transcriptoma y asignan lecturas a transcritos de forma probabilística. Ventajas: mucha velocidad y menor uso de recursos; limitaciones: menos información posicional que un mapeo tradicional, pero suficiente para cuantificación y muchos análisis de expresión.",
  },
  {
    pattern:
      "^(que es (el )?bias gc|que significa bias gc|que es sesgo gc|explicame bias gc|como afecta (el )?contenido gc a la secuenciacion)[!., ]*$",
    response:
      "⚖️ El bias GC es un sesgo técnico donde regiones con alto o bajo contenido en GC se representan de forma desigual en las lecturas. Puede afectar la eficiencia de PCR, la captura de librerías y la secuenciación, provocando que ciertos genes parezcan más o menos expresados por razones técnicas. Se detecta en controles de calidad y se puede corregir con métodos estadísticos o usando herramientas que modelan el sesgo en la cuantificación.",
  },
  {
    pattern:
      "^(que es (la )?normalizacion por (tamano|tamanio) de libreria|explicame normalizacion por (tamano|tamanio) de libreria|para que sirve normalizacion por (tamano|tamanio) de libreria)[!., ]*$",
    response:
      "⚖️ La normalización por tamaño de librería ajusta los conteos para tener en cuenta diferencias en el número total de lecturas entre muestras. Su objetivo es que las comparaciones de expresión reflejen cambios biológicos y no diferencias en profundidad de secuenciación. Métodos simples dividen por el total de lecturas (o por millones de lecturas) y métodos más sofisticados (TMM, RLE) corrigen además composiciones sesgadas.",
  },
  {
    pattern:
      "^(que es (la )?dispersion en rnaseq|que significa dispersion en rnaseq|explicame dispersion en rnaseq|para que sirve (la )?dispersion en rnaseq)[!., ]*$",
    response:
      "📈 La dispersión en RNA-Seq mide la variabilidad biológica y técnica entre replicados para un mismo gen, normalmente en el contexto de modelos de conteo (por ejemplo, distribución binomial negativa). Estimar correctamente la dispersión es clave para evaluar la significancia de cambios de expresión: dispersión alta reduce la potencia estadística y dispersión baja aumenta la confianza en diferencias observadas.",
  },
  {
    pattern:
      "^(que es (el )?ajuste por multiples pruebas|que es correccion por comparaciones multiples|explicame ajuste por multiples pruebas|para que sirve ajustar por multiples pruebas)[!., ]*$",
    response:
      "🧾 El ajuste por múltiples pruebas corrige el problema de falsos positivos cuando se realizan miles de tests simultáneos (por ejemplo, un test por cada gen). Métodos como Bonferroni o FDR (Benjamini-Hochberg) controlan la tasa de error para que las conclusiones sean más robustas. Sin esta corrección, muchas diferencias aparentemente significativas serían en realidad debidas al azar.",
  },
  {
    pattern:
      "^(que es (un )?pvalor ajustado|que significa pvalor ajustado|que es adjusted pvalue|explicame pvalor ajustado|para que sirve (el )?pvalor ajustado)[!., ]*$",
    response:
      "✅ El pvalor ajustado es el resultado del pvalor original tras aplicar una corrección por múltiples pruebas (por ejemplo, FDR). Refleja la probabilidad de obtener el resultado bajo la hipótesis nula considerando el número total de comparaciones. Usar pvalores ajustados reduce falsos positivos y ayuda a priorizar genes con evidencia estadística real.",
  },
  {
    pattern:
      "^(que es (el )?fc|que es (el )?fold change|que significa (el )?fold change|que es cambio de plegamiento|explicame fold change|para que sirve fold change)[!., ]*$",
    response:
      "🔁 El fold change (o ratio) mide el cambio relativo en expresión entre dos condiciones (por ejemplo, condición A vs B). Suele expresarse como log2 fold change para simetrizar aumentos y disminuciones. Complementado con la significancia estadística, permite identificar genes con cambios biológicamente relevantes.",
  },
  {
    pattern:
      "^(que es (el )?clustering jerarquico|que significa (el )?clustering jerarquico|explicame clustering jerarquico|para que sirve clustering jerarquico|que es (el )?agrupamiento jerarquico)[!., ]*$",
    response:
      "🌳 El clustering jerárquico es un método para agrupar genes o muestras según su similitud de expresión, construyendo una jerarquía de clústeres mediante enlaces sucesivos. En heatmaps se usa para ordenar filas y columnas, facilitando la identificación de grupos coexpresados y patrones de respuesta. Elegir la medida de distancia y el método de enlace influye en el resultado, por eso es importante justificar la elección.",
  },
  {
    pattern:
      "^(que es (un )?pipeline (de )?rnaseq|que significa (un |el )?pipeline( de)? rnaseq|que es un flujo de trabajo rnaseq|explicame pipeline rnaseq|para que sirve (un )?pipeline rnaseq)[!., ]*$",
    response:
      "🛠 Un pipeline RNA-Seq es la secuencia organizada de pasos que transforma datos crudos (FASTQ) en resultados biológicos interpretables: control de calidad, trimming, alineamiento o cuantificación, cuantificación de genes/transcritos, normalización, análisis diferencial y anotación funcional. Automatizar el pipeline mejora la reproducibilidad, trazabilidad y facilita la comparación entre experimentos.",
  },
];
