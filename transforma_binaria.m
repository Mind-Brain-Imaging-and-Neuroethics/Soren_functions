function s=transforma_binaria(serie,num_simbolos,umbral)

% Funcion:   transformar_senal
%
%                   Transforma una se�al temporal en una serie compuesta por un
%                   n�mero de s�mbolos determinado
%
% Argumentos de entrada:
%   serie:          se�al temporal a transformar
%   num_simbolos:   n�mero de s�mbolos a utilizar en la transformaci�n (2 � 3)
%   umbral:         umbral de decisi�n a aplicar en la transformaci�n
% Argumentos de salida:
%   s:              se�al transformada
%
% Ultima modificacion:   11 de Diciembre de 2007
%
% Autor:            Alicia Rodrigo de Diego y Jose Victor Marcos Martin
% Modificado por:   Daniel �lvarez Gonz�lez
%

    if num_simbolos==2
        if strcmp(umbral,'mediana')  
            mediana=median(serie);
            for i=1:1:length(serie)
                if serie(i)<mediana
                    s(i)=0;
                else
                    s(i)=1;
                end
            end
        else
            media=mean(serie);
            for i=1:1:length(serie)
                if serie(i)<media
                    s(i)=0;
                else
                    s(i)=1;
                end
            end
        end
        
    else
        mediana=median(serie);
        mx=abs(max(serie));
        mn=abs(min(serie));
        
        td1=mediana-mn/16;
        td2=mediana+mx/16;
        
        for i=1:1:length(serie)
            if serie(i)<=td1
                s(i)=0;
            else if serie(i)<td2
                    s(i)=1;
                else
                    s(i)=2;
                end
            end
        end
    end