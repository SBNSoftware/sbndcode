echo "\documentclass{article}
\usepackage[a4paper, margin=3cm]{geometry}
\usepackage{graphicx}
\usepackage{pgffor}
\usepackage[hidelinks]{hyperref}

\title{SBND CRT Channel Mapping Displays \\\\ \vspace{1em} \small \textit{Produced using} \texttt{sbndcode v09\_82\_02\_01} \textit{\&} \texttt{sbnd\_v02\_00.gdml}}
\author{Henry Lay \\\\ \small h.lay@lancaster.ac.uk}" > crt_channel_mapping_evds.tex

walls=(bottom south north west east toplow tophigh)
wallnames=(Bottom South North West East "Top Low" "Top High")

for wall in "${walls[@]}"
do
    list=$(ls /exp/sbnd/data/users/hlay/crt_channel_mapping/${wall}_wall/*_front.pdf)
    echo -n "\newcommand*{\\"$wall"ids}{" >> crt_channel_mapping_evds.tex

    for item in ${list}
    do
	name=$(echo $item | cut -d '/' -f 9)
	number=$(echo $name | cut -d '_' -f 2)
	echo -n $number, >> crt_channel_mapping_evds.tex
    done

    sed -i '$ s/.$//' crt_channel_mapping_evds.tex
    echo "}" >> crt_channel_mapping_evds.tex
done

echo "\begin{document}

\maketitle

\centering
\vspace{2em}

\includegraphics[width=.6\textwidth]{/exp/sbnd/data/users/hlay/crt_channel_mapping/luphysics_logo.png}

\vspace{2em}

\includegraphics[width=.5\textwidth]{/exp/sbnd/data/users/hlay/crt_channel_mapping/sbnd_pride_transparent.png}
\flushleft
\newpage
\tableofcontents
\newpage
\section{Explanation}
This document contains a series of illustrations created using the \texttt{CRTEventDisplay} tool originally written by Tom Brooks \& heavily developed by myself. It shows the position of the various CRT modules according to the gdml file used in SBND simulation and reconstruction. The document is split into sections for the different tagger walls. For each module three illustrations are provided: front, top and side views. The axes show detector coordinates (X, Y and Z) and \`\`building coordinates\" (South, West and Up). The relevant module is shown in green. The TPCs are shown in grey in the centre for reference. The black outer is the full tagger wall. The thin grey are other modules in the wall. The red is the FEB position and the blue corresponds to the end of the FEB with channel 0 (the ethernet ports).
" >> crt_channel_mapping_evds.tex

for i in "${!walls[@]}"
do
    echo "\newpage
\section{${wallnames[i]} Wall}
\begingroup
\foreach \x in \\${walls[i]}ids
{
    \newpage
    \subsection{volCRTModule\x\_\x}
    \begin{center}
            \includegraphics[width=.85\textwidth]{/exp/sbnd/data/users/hlay/crt_channel_mapping/${walls[i]}_wall/volCRTModule\x_\x_front.pdf}\\\\
	    \includegraphics[width=.85\textwidth]{/exp/sbnd/data/users/hlay/crt_channel_mapping/${walls[i]}_wall/volCRTModule\x_\x_top.pdf}\\\\
	    \includegraphics[width=.85\textwidth]{/exp/sbnd/data/users/hlay/crt_channel_mapping/${walls[i]}_wall/volCRTModule\x_\x_side.pdf}
    \end{center}
}
\endgroup" >> crt_channel_mapping_evds.tex
done

echo "\end{document}" >> crt_channel_mapping_evds.tex

pdflatex --shell-escape -output-directory /exp/sbnd/data/users/hlay/crt_channel_mapping/tex_work crt_channel_mapping_evds.tex
pdflatex --shell-escape -output-directory /exp/sbnd/data/users/hlay/crt_channel_mapping/tex_work crt_channel_mapping_evds.tex
