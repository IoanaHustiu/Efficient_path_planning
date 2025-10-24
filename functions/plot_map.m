function h = plot_map(B, varargin)
%PLOT_MAP  Muestra un mapa binario en coordenadas cartesianas (0,0 abajo-izquierda).
%
%   h = PLOT_MAP(B) visualiza la matriz binaria B (H×W), donde 1=obstáculo y 0=libre,
%   usando un sistema de coordenadas cartesiano: el origen (0,0) está en la esquina
%   inferior izquierda, x crece a la derecha e y hacia arriba.
%
%   B debe estar ya en convención cartesiana si proviene de un *loader* que
%   hace flipud (por ejemplo, tu `load_map_and_scen`). Si NO es el caso, puedes
%   pasar 'BIsCartesian', false y la función aplicará flipud para visualizarlo.
%
%   h = PLOT_MAP(B, Name,Value,...) permite ajustar opciones:
%       'BIsCartesian'  - (logical)  true por defecto. Si false, se aplica flipud(B).
%       'ShowGrid'      - (logical)  true para mostrar cuadrícula (default=false).
%       'GridAlpha'     - (double)   transparencia de la cuadrícula [0..1] (default=0.07).
%       'Axes'          - (axes)     ejes existentes donde dibujar (default=[] crea figura).
%       'Title'         - (char/str) título de la figura (default: descriptivo).
%       'ShowColorbar'  - (logical)  muestra barra de color (default=false).
%
%   Salida:
%       h.fig  - handle de la figura (si se creó).
%       h.ax   - handle de los ejes.
%       h.img  - handle del objeto image (imagesc).
%
%   Ejemplos:
%       % B ya en coordenadas cartesianas (p. ej., tras load_map_and_scen):
%       h = plot_map(B, 'ShowGrid', true);
%
%       % Si B aún está en convención "imagen" (0,0 arriba-izquierda):
%       h = plot_map(B, 'BIsCartesian', false, 'ShowGrid', true);
%
%   Notas:
%     - Se representa ~B para que las celdas libres aparezcan en blanco y los
%       obstáculos en negro con colormap(gray).
%     - El eje y se fuerza con 'YDir','normal' para que crezca hacia arriba.
%
%   Autor: Cristian Mahulea | Fecha: October 20, 2025

% --------------------
% Validación de entrada
% --------------------
validateattributes(B, {'numeric','logical'}, {'2d','nonempty'}, mfilename, 'B', 1);
if ~islogical(B)
    % Tolerar 0/1 numéricos
    B = B ~= 0;
end

% Parámetros opcionales
ip = inputParser;
ip.addParameter('BIsCartesian', true,  @(x)islogical(x)&&isscalar(x));
ip.addParameter('ShowGrid',     false, @(x)islogical(x)&&isscalar(x));
ip.addParameter('GridAlpha',    0.07,  @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
ip.addParameter('Axes',         [],    @(x)isempty(x) || isa(x,'matlab.graphics.axis.Axes'));
ip.addParameter('Title',        [],    @(x)ischar(x)||isstring(x)||isempty(x));
ip.addParameter('ShowColorbar', false, @(x)islogical(x)&&isscalar(x));
ip.parse(varargin{:});

BIsCartesian = ip.Results.BIsCartesian;
showGrid     = ip.Results.ShowGrid;
gridAlpha    = ip.Results.GridAlpha;
ax           = ip.Results.Axes;
ttl          = ip.Results.Title;
showCbar     = ip.Results.ShowColorbar;

% Convertir a cartesiano si hace falta
if ~BIsCartesian
    B = flipud(B);
end

[H, W] = size(B);

% --------------------
% Crear figura / ejes
% --------------------
createdFig = false;
if isempty(ax) || ~isvalid(ax)
    h.fig = figure('Name','Map (Cartesian coordinates)','Color','w');
    ax = axes('Parent', h.fig);
    createdFig = true;
else
    h.fig = ancestor(ax, 'figure');
end
h.ax = ax;

% --------------------
% Dibujo
% --------------------
axes(ax); %#ok<LAXES>
cla(ax); hold(ax,'on');

% Mostramos libre=blanco, obstáculo=negro -> usar ~B con colormap(gray)
h.img = imagesc(ax, [0 W-1], [0 H-1], ~B);
colormap(ax, gray);
set(ax, 'YDir','normal');     % y hacia arriba
axis(ax, 'equal'); axis(ax, 'tight'); box(ax, 'on');
xlabel(ax, 'x'); ylabel(ax, 'y');

% Título
if isempty(ttl)
    ttl = 'Black = obstacle, White = free';
end
title(ax, ttl, 'Interpreter','none');

% Cuadrícula opcional (en líneas de celda)
if showGrid
    set(ax, 'XTick', 0:W, 'YTick', 0:H, 'GridAlpha', gridAlpha);
    grid(ax, 'on');
else
    % ticks más livianos
    set(ax, 'XTickMode','auto', 'YTickMode','auto');
end

% Colorbar opcional
if showCbar
    colorbar(ax);
end

% Traer al frente cualquier overlay futuro
uistack(h.img,'bottom');

% Si se creó figura, optimizar renderizado
if createdFig
    drawnow;
end
end
