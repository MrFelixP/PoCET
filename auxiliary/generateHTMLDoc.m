function generateHTMLDoc(filename, output_directory, title)
% generateHTMLDoc(filename, [output_directory, [file_extension, title]])
%
% generateHTMLDOC Generates a HTML documentation out of the help text
% of an m-File.
%
% filename: Name of the m-File
% output_directory (optional, default: '.'): Output directory of the html-file
% title (optional): Title of the HTML-File. If not given, then the title is
% determined by the name of the m-file. The title can also be set in the
% documentation comment using a double percent sign.
% Example:
% %% This is the title
% % This is the other text.
%
% See the comments in the source code for details.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ========================
% HOW TO USE THIS FUNCTION
% ========================
% Write the documentation in form of a comment in an m-file. The first
% comment found in the m-file will be considered for the generation of the
% html file.
% Some things will be formatted automatically, otherwise you can use html
% code to format the text, include figures, links or lists.
% Examples for such html formatting options are
%   This is<b>bold text</b>.
%   This is<i>italic text</i>.
%   <img src="path to image" /> includes an image
%   This creates a <a href="target of the link">link</a> to another file.
%   This creates a list with bullet points:
%   <ul>
%     <il>list item 1</il>
%     <il>list item 2</il>
%     <il>list item 3</il>
%   </ul>
%
%
%
% The following things will be formatted automatically by generateHTMLDoc.
% 
% 1. Leading percent signs of the comment will be deleted.
% 2. Text in the form
%    =======
%    Section
%    =======
%    and
%    ----------
%    Subsection
%    ----------
%    will be formatted as big section titles or smaller subsection titles.
% 3. EVERY word starting with PoCET (case sensitive, e.g. PoCETcompose,
%    PoCETestimate but also PoCETfoobar1234_blubb) will be replaced by a
%    link to the corresponding html-file (e.g. PoCETcompose <a
%    href="PoCETcompose.html">PoCETcompose</a>).
%    This has the following side effect: If you manually link to a html
%    file like PoCETToolbox_product_page.html, then put the href in
%    double quotes. This means, use the form
%    <a href="PoCETToolbox_product_page.html">Link text</a> instead of
%    <a href=PoCETToolbox_product_page.html>Link text</a>. In the latter
%    case the string PoCETToolbox_product_page will be replaced by
%    <a href="PoCETToolbox_product_page.html">PoCETToolbox_product_page</a>
%    resulting in invalid html code. The double quotes prevent this
%    behaviour.
% 4. PoCEToption-names like ESTIMATE.bisectioning.use will be automatically
%    replaced by a link to the according documentation page for the
%    specific option. E.g. ESTIMATE.bisectioning.use becomes
%    <a href="options_ESTIMATE.html#ESTIMATE.bisectioning.use">ESTIMATE.bisectioning.use</a>
% 5. The documentation of PoCEToptions is handled in a special way. The
%    help text of PoCEToptions also contains the documentation for all
%    options. generateHTMLDoc splits the documentation in several files,
%    such that a file PoCEToptions.html containing the documentation of the
%    PoCEToptions class and several files like options_COMPOSE.html or
%    options_DISPLAY.html containing the documentation for the specific
%    options will be created.
% 6. Greater than and lower than signs will be replaced. Otherwise invalid
%    html code could be generated. The following replacements are carried
%    out (excluding the single quotes but including whitespaces):
%      ' < '  => ' &lt; ' 
%      '<='   => '&lt;='
%      ' > '  => ' &gt; '
%      ' >= ' => ' &gt;= '
% 7. Mark ---- EXPERIMENTAL ---- in bold and red.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(filename(end-1:end), '.m')
    warning('PoCET:http_generation:no_m_file', '%s is not an m-file', filename);
end

if nargin < 2
    output_directory = '.';
end

file_extension = '.html';

filecontent = regexprep(fileread(filename), '\r\n', '\n');
% get function name or class name via regular expression
% name = regexp(filecontent2, '^(?:(function (?:\[?[a-zA-Z0-9_\s,]*?\]? *= *)?)|(classdef ))([a-zA-Z0-9]+)', 'tokens');                                  
% if isempty(name)
%     return
% end

% get function name or class name via filename
name = regexp(filename, '(?:.*/)*(.+)\.m', 'tokens');
name = name{1}{1};

if nargin < 3
    if length(filecontent) >= 2 && strcmp(filecontent(1:2), '%%')
        title = regexp(filecontent, '%% ?(.+?)\n', 'tokens');
        title = title{1}{1};
        filecontent = regexp(filecontent, '.*?\n(.*)', 'tokens');
        filecontent = filecontent{1}{1};
    else
        title = name;
    end
end

% extract docstring from filecontent
docstring = regexp(filecontent, '(% ?.+?\n)+', 'tokens');

if isempty(docstring)
    % stop if the file does not contain a docstring
    warning('PoCET:http_generation:no_docstring', '%s does not contain a documentation comment. No html file has been created.', filename);
    return;
end

% remove leading percent signs and whitespace from docstring
docstring = regexprep(docstring{1}{1}, '(\n% ?)|(^% ?)', '\n');
if ~strcmp(name, 'PoCEToptions')
    html_doc = create_html(docstring, title);
    html_file = fopen(fullfile(output_directory, [name file_extension]), 'w');
    fwrite(html_file, html_doc);
    fclose(html_file);
else
    createOptionsDocumentation(docstring, output_directory, file_extension);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function html_source = create_html(text, title)
% create_html
% Creates the HTML source out of the help text using a template file.
html_source = ['<html>' ...
'<head>' ...
'	<style type="text/css">' ...
'		h1 {color: #900; font-size: 24px; padding: 2px 0px 2px 0px; margin: 10px 0px 0px 5px; vertical-align: middle;}' ...
'		h2 {color: #900; font-weight: bolder}' ...
'		.doctext  {display: block; font-family: monospace; white-space: pre; margin: 1em 0px;}' ...
'		.keyword  {color: #900; font-weight: bold;}' ...
'		.header {background-color: #900; color: white; font-weight: bold}' ...
'	</style>' ...
'</head>' ...
'<title>' ...
'<!-- TITLE PLACEHOLDER -->' ...
'</title>' ...
'<body>' ...
'<div class="header">PoCET documentation: <!-- TITLE PLACEHOLDER --></div>' ...
'<h1><!-- TITLE PLACEHOLDER --></h1>' ...
'<pre class="doctext"><!-- TEXT PLACEHOLDER --></pre>' ...
'</body>' ...
'</html>'];
% html_source = fileread('__template.html');
html_source = strrep(html_source, '<!-- TITLE PLACEHOLDER -->', title);

% process help text

% mask lower than and greater than signs, the Matlab browser misunderstands
% as tags
text = strrep(text, ' < ', ' &lt; '); 
text = strrep(text, '<=', '&lt;=');
text = strrep(text, ' > ', ' &gt; ');
text = strrep(text, ' >= ', ' &gt;= ');

% search for 
%  =========
%  SECTIONS:
%  =========
[sections, sectionnames] = regexp(text, '\n==+\n(.+?):?\n==+', 'match', 'tokens');
for i = 1:length(sections)
    text = strrep(text, sections{i}, ['</pre><h2>' sectionnames{i}{1} '</h2><pre class="doctext">']);
end

% search for
%  ------------
%  SUBSECTIONS:
%  ------------
[subsections, subsectionnames] = regexp(text, '\n--+\n(.+?):?\n--+', 'match', 'tokens');
for i = 1:length(subsections)
    text = strrep(text, subsections{i}, ['</pre><h3>' subsectionnames{i}{1} '</h3><pre class="doctext">']);
end

% make Inputs: and Returns: bold
text = strrep(text, 'Inputs:', '<b>Inputs</b>');
text = strrep(text, 'Returns:', '<b>Returns</b>');

% create links to other PoCET-Functions
tokens = regexp(text, '([A-Za-z0-9_]*PoCET[A-Za-z0-9_]*)', 'tokens');
admit_keywords = {};
for i = 1:length(tokens)
    admit_keywords{end+1} = tokens{i}{1}; %#ok<AGROW>
end
admit_keywords = unique(admit_keywords);
admit_keywords = setdiff(admit_keywords,'PoCET');

% also create automatic links to installPoCET and displayMessage
admit_keywords{end+1} = 'PoCETsimMonteCarlo';
admit_keywords{end+1} = 'PoCETwriteFiles';

admit_keywords = sort_by_decreasing_length(admit_keywords);

for i = 1:length(admit_keywords)
    if strcmp(admit_keywords{i}, title)
        continue
    end
    text = mystrrep(text, admit_keywords{i}, sprintf('<a href="%s.html">%s</a>', admit_keywords{i}, admit_keywords{i}));
end

% % create links to PoCET options when we don't generate a options
% % documentation page
% if isempty(strfind(title, 'PoCEToptions: '))
%     optionsprefixes = {'YALMIP', 'COMPOSE', 'ESTIMATE', 'PLOT', 'PARALLEL', 'DISPLAY', 'SIMULATE'};
%     for idx1 = 1:length(optionsprefixes)
%         %     tokens = regexp(text, sprintf('(%s\..+)\W', optionsprefixes{idx1}), 'tokens');
%         tokens = regexp(text, ['(' optionsprefixes{idx1}, '[', regexptranslate('escape', '.'), 'A-Za-z]+?)[''\s]'], 'tokens');
%         options_names = {};
%         for idx2 = 1:length(tokens)
%             options_names{end+1} = tokens{idx2}{1}; %#ok<AGROW>
%         end
%         options_names = unique(options_names);
%         options_names = sort_by_decreasing_length(options_names);
%         for idx2 = 1:length(options_names)
%             token = options_names{idx2};
%             text = mystrrep(text, token, sprintf('<a href="options_%s.html#%s">%s</a>', optionsprefixes{idx1}, token, token));
%         end
%     end
% end
% % mark the function or class name
% % text = regexprep(text, ['(' title ')\W'], sprintf('<b class="keyword">%s</b>', title), 'tokens');
% % text = strrep(text, title, ['<b class="keyword">' title '</b>']);
% 
% % Mark ---- EXPERIMENTAL ---- in red and bold
% text = strrep(text, '---- EXPERIMENTAL ----', '<b class="header">---- EXPERIMENTAL ----</b>');

html_source = strrep(html_source, '<!-- TEXT PLACEHOLDER -->', text);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = sort_by_decreasing_length(a)
r = a;
% very simple bubblesort
for idx1 = 1:length(r)
    for idx2 = 1:length(r)-1
        len1 = length(r{idx2});
        len2 = length(r{idx2+1});
        if len1 < len2
            tmp = r{idx2};
            r{idx2} = r{idx2+1};
            r{idx2+1} = tmp;
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modifiedstr = mystrrep(origstr, oldsubstr, newsubstr)
% This works like MATLAB's strrep function, but it won't replace the substr
% if the character before the substring ist a > or " or if the character 
% after the substring ist alphanumeric. This is done to avoid the
% replacement of strings which were already replaced by html-a-tags.
% Example: Assume ESTIMATE.outerBounding.useOptimization was already
% replaced by the automatic link <a
% href="...">ESTIMATE.outerBounding.useOptimization</a>/ Then a simple
% newtext = strrep(text, 'ESTIMATE.outerBounding.use', ...) would break
% this link. mystrrep avoids this.
oldsubstr_idx = strfind(origstr, oldsubstr);
modifiedstr = origstr;

% test if origstr contains oldsubstr
if ~isempty(oldsubstr_idx)
    % position of the first oldsubstr in origstr
    oldsubstr_idx = oldsubstr_idx(1);
    
    % split origstr before and after oldsubstr
    presub = origstr(1:oldsubstr_idx-1);
    postsub = origstr(oldsubstr_idx + length(oldsubstr):length(origstr));
    
    % test if replacement is allowed
    if (oldsubstr_idx > 1 && ~isempty(strfind('>"', presub(end)))) || ...
        (~isempty(postsub) && isalpha_num(postsub(1)))
        % it's not allowed
        modifiedstr = [presub, oldsubstr, mystrrep(postsub, oldsubstr, newsubstr)];
    else
        % it's allowed
        modifiedstr = [presub, newsubstr, mystrrep(postsub, oldsubstr, newsubstr)];
    end
end
end

function r = isalpha_num(s)
r = ~isempty(strfind('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ', upper(s)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function createOptionsDocumentation(docstring, output_directory, file_extension)
% create the documentation for PoCEToptions
PoCEToptions_doctext = regexp(docstring, '(.*?)\n==+', 'tokens');
html_doc = create_html(PoCEToptions_doctext{1}{1}, 'PoCEToptions');
html_file = fopen(fullfile(output_directory, ['PoCEToptions' file_extension]), 'w');
fwrite(html_file, html_doc);
fclose(html_file);

% create the documentations for the options categories
prefixes = regexp(docstring, '\n==+\nPrefix (.+?):?\n==+', 'tokens');
for prefixidx = 1:length(prefixes)
    prefix = prefixes{prefixidx}{1};
    title = sprintf('PoCEToptions: %s', prefix);
    if prefixidx == length(prefixes)
        prefixdocstring = regexp(docstring, sprintf('\n==+\nPrefix %s:?\n==+\n(.*)\n', prefix), 'tokens');
    else
        prefixdocstring = regexp(docstring, sprintf('\n==+\nPrefix %s:?\n==+\n(.*?)\n==+', prefix), 'tokens');
    end
    prefixdocstring = prefixdocstring{1}{1};
    % bold option names
    optionnames = regexp(prefixdocstring, ['^(', prefix, '\..*?) '], 'tokens', 'lineanchors');
    optionnames = sort_by_decreasing_length(optionnames);
    for optionidx = 1:length(optionnames)
        prefixdocstring = mystrrep(prefixdocstring, optionnames{optionidx}{1}, sprintf('<b id="%s">%s</b>', optionnames{optionidx}{1}, optionnames{optionidx}{1}));
    end
    html_doc = create_html(prefixdocstring, title);
    html_file = fopen(fullfile(output_directory, ['options_' prefix file_extension]), 'w');
    fwrite(html_file, html_doc);
    fclose(html_file);
end

% create the overview page
template = 'The following option categories are available:\n\n%s';

list_string = '';
for i = 1:length(prefixes)
    prefix = prefixes{i}{1};
    list_string = sprintf('%s<li><a href=options_%s%s>%s</a></li>', list_string, prefix, file_extension, prefix);
end
list_string = sprintf('<ul>%s</ul>', list_string);
html_source = create_html(sprintf(template, list_string), 'Options reference');
html_file = fopen(fullfile(output_directory, ['options_overview' file_extension]), 'w');
fwrite(html_file, html_source);
fclose(html_file);
end

