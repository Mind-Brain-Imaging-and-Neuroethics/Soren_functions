function [LZC,meanLZC]=LZ_Complexity_MultiBin_Window(registro,tam,solapamiento)
%   function LZcomplexity=lzcomplexity_tramas(registro,umbral_db,tam,solapamiento)
%
%   La serie de entrada es dividida en segmentos, de forma que se obtiene 
%   un valor de LZC para cada segmento. Finalmente se realiza un 
%   promediado para obtener un único valor de LZC para la serie.
%%
%   Argumentos:
%
%   registro:       Serie de datos de entrada de la que estimaremos su LZC.
%   umbral_db:      Umbral empleado en la conversión binaria de la serie.
%                   Tomará los valores 'media' o 'mediana'.
%   num_simbolos:   Número de símbolos empleados en la conversión binaria
%                   de la serie de datos original.
%   tam:            Número de muestras de los segmentos en que se dividirá la 
%                   serie de datos original. Si tam = 0, se aplicará el
%                   algoritmo sobre la serie original sin segmentar.
%   solapamiento:   valor entre 0 y 1 que indica el tanto por uno de solapamiento
%
%   Variables de salida:
%   LZC:        Valor de LZC de la serie de datos de entrada.
%
%
%   Última actualización: 17/01/2013
%   Autor: Carlos Gómez Peña
%

% Calculamos el número de veces a aplicar el algoritmo según el número de 
% segmentos de longitud tam muestras en los que se puede dividir la serie 
% de datos de entrada.

warning off

if tam==0   % Si tam = 0 el algoritmo se calcula sobre el registro completo
    n_tramas=1;
    tam=length(registro);
    solapamiento=0; % Independientemente de lo que haya introducido el usuario
else
    n_tramas=fix((length(registro)-tam)/(tam*(1-solapamiento))+1);
end

% Calculamos LZC para cada trama
for l=1:n_tramas
    % Extraemos la trama correspondiente del registro original.
    if l==n_tramas
        trama=registro((1+(l-1)*tam*(1-solapamiento)):end);
    else
        trama=registro((1+(l-1)*tam*(1-solapamiento)):((l-1)*tam*(1-solapamiento)+tam));
    end
    
    % Número de muestras del segmento.
    n=length(trama);
    
    LZC(l) = LZC_complexity_V2(trama);
    
    meanLZC = mean(LZC);
end

