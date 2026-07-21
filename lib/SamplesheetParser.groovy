import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

class SamplesheetParser {

    static Map parse(def inputValue, def comparisonValue) {
        Path sheet = Paths.get(inputValue.toString()).toAbsolutePath().normalize()
        if (!Files.exists(sheet) || !Files.isRegularFile(sheet)) {
            throw new IllegalArgumentException("Samplesheet does not exist or is not a regular file: ${sheet}")
        }

        List<String> lines = Files.readAllLines(sheet)
            .findAll { it != null && !it.trim().isEmpty() && !it.trim().startsWith('#') }
        if (lines.size() < 2) {
            throw new IllegalArgumentException("Samplesheet must contain a header and at least one data row: ${sheet}")
        }

        char delimiter = detectDelimiter(lines[0])
        List<String> headers = splitDelimited(lines[0], delimiter)
            .collect { it.replace('\uFEFF', '').trim() }

        Map<String,Integer> headerIndex = [:]
        headers.eachWithIndex { h, i -> headerIndex[h.toLowerCase()] = i }
        Map<String,List<String>> aliases = [
            id        : ['samples', 'sample', 'sample_id'],
            r1        : ['r1'],
            r2        : ['r2'],
            condition : ['condition'],
            replicate : ['replicate']
        ]
        Map<String,Integer> columns = [:]
        aliases.each { canonical, names ->
            def found = names.find { headerIndex.containsKey(it) }
            if (found == null) {
                throw new IllegalArgumentException(
                    "Missing required samplesheet column for '${canonical}'. " +
                    "Accepted header(s): ${names.join(', ')}. Observed: ${headers.join(', ')}"
                )
            }
            columns[canonical] = headerIndex[found]
        }

        Path baseDir = sheet.parent ?: Paths.get('.').toAbsolutePath().normalize()
        List<Map> records = []
        Set<String> seenIds = [] as Set

        lines.drop(1).eachWithIndex { line, idx ->
            int lineNumber = idx + 2
            List<String> fields = splitDelimited(line, delimiter)
            if (fields.size() != headers.size()) {
                throw new IllegalArgumentException(
                    "Samplesheet line ${lineNumber} has ${fields.size()} field(s), but header has ${headers.size()}: ${line}"
                )
            }
            Map rec = [
                id        : fields[columns.id].trim(),
                r1        : fields[columns.r1].trim(),
                r2        : fields[columns.r2].trim(),
                condition : fields[columns.condition].trim(),
                replicate : fields[columns.replicate].trim()
            ]

            ['id','r1','r2','condition','replicate'].each { key ->
                if (!rec[key]) {
                    throw new IllegalArgumentException("Empty '${key}' value on samplesheet line ${lineNumber}")
                }
            }
            validateToken(rec.id, 'Samples', lineNumber)
            validateToken(rec.condition, 'Condition', lineNumber)
            validateToken(rec.replicate, 'Replicate', lineNumber)

            if (!seenIds.add(rec.id)) {
                throw new IllegalArgumentException("Duplicate Samples value '${rec.id}' on line ${lineNumber}")
            }

            rec.r1 = resolveInputPath(baseDir, rec.r1, lineNumber, 'R1').toString()
            rec.r2 = resolveInputPath(baseDir, rec.r2, lineNumber, 'R2').toString()
            if (rec.r1 == rec.r2) {
                throw new IllegalArgumentException("R1 and R2 resolve to the same file for sample '${rec.id}'")
            }
            records << rec
        }

        String comparison = comparisonValue.toString().trim()
        String[] contrast = comparison.split(':', -1)
        if (contrast.size() != 2 || !contrast[0].trim() || !contrast[1].trim()) {
            throw new IllegalArgumentException(
                "--danpos_comparison must contain exactly two non-empty condition names separated by one colon, e.g. 'iKO_24h:noKO_24h'"
            )
        }
        String treatment = contrast[0].trim()
        String control = contrast[1].trim()
        validateToken(treatment, 'DANPOS treatment condition', 0)
        validateToken(control, 'DANPOS control condition', 0)
        if (treatment == control) {
            throw new IllegalArgumentException("DANPOS treatment and control conditions must be different")
        }

        Set<String> conditions = records.collect { it.condition } as Set
        [treatment, control].each { condition ->
            if (!conditions.contains(condition)) {
                throw new IllegalArgumentException(
                    "Condition '${condition}' from --danpos_comparison is absent from the Condition column. " +
                    "Available conditions: ${conditions.sort().join(', ')}"
                )
            }
        }

        return [
            records    : records,
            treatment  : treatment,
            control    : control,
            conditions : conditions.sort(),
            delimiter  : delimiter.toString()
        ]
    }

    private static Path resolveInputPath(Path baseDir, String rawPath, int lineNumber, String column) {
        String expanded = rawPath
        if (expanded.startsWith('~/')) {
            expanded = System.getProperty('user.home') + expanded.substring(1)
        }
        Path p = Paths.get(expanded)
        if (!p.isAbsolute()) p = baseDir.resolve(p)
        p = p.toAbsolutePath().normalize()
        if (!Files.exists(p) || !Files.isRegularFile(p)) {
            throw new IllegalArgumentException("${column} file on samplesheet line ${lineNumber} does not exist: ${p}")
        }
        return p
    }

    private static void validateToken(String value, String label, int lineNumber) {
        if (!(value ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
            String where = lineNumber > 0 ? " on line ${lineNumber}" : ''
            throw new IllegalArgumentException(
                "${label}${where} contains unsupported characters: '${value}'. " +
                "Use letters, numbers, '.', '_' or '-', beginning with a letter or number."
            )
        }
    }

    private static char detectDelimiter(String header) {
        int tabs = header.count('\t')
        int commas = header.count(',')
        if (tabs == 0 && commas == 0) {
            throw new IllegalArgumentException("Could not detect CSV/TSV delimiter from header: ${header}")
        }
        return tabs >= commas ? ('\t' as char) : (',' as char)
    }

    private static List<String> splitDelimited(String line, char delimiter) {
        List<String> out = []
        StringBuilder current = new StringBuilder()
        boolean quoted = false
        for (int i = 0; i < line.length(); i++) {
            char c = line.charAt(i)
            if (c == '"' as char) {
                if (quoted && i + 1 < line.length() && line.charAt(i + 1) == ('"' as char)) {
                    current.append('"')
                    i++
                } else {
                    quoted = !quoted
                }
            } else if (c == delimiter && !quoted) {
                out << current.toString()
                current.setLength(0)
            } else {
                current.append(c)
            }
        }
        if (quoted) {
            throw new IllegalArgumentException("Unclosed quoted field in samplesheet line: ${line}")
        }
        out << current.toString()
        return out
    }
}
