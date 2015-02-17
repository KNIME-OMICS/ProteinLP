package uni.tubingen.protein;

import java.io.PrintWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;

/**
 * Date: Sep 14, 2007
 * Time: 21:52:04 PM
 *
 */
public class Protein
{
    private String m_header;  // the header information of the uni.tubingen.protein is stored here
    private byte[] m_bytes;   // the uni.tubingen.protein sequence is stored here
    private double m_mass;
    private String m_lookup;  // the unique ID of uni.tubingen.protein in the fasta DB

    private String m_origHeader;  // original line in DB for head information of uni.tubingen.protein 
    private Map m_identifierMap;

    //known identifier types.  Multiple identifiers found in fasta files can often
    //boil down to the same thing
    public static HashMap IdentTypeMap = new HashMap();

    /* for parsing header lines of FASTA files */
    public static final String SEPARATOR_PATTERN = "\\|";
    public static final String SEPARATOR_CHAR = "|";

    //populate the hashmap of known identifier types
    static
    {
        IdentTypeMap.put("GI", "GI");
        IdentTypeMap.put("REF", "REFSEQ");
        IdentTypeMap.put("GB", "Genbank");
        IdentTypeMap.put("EMB", "Genbank");
        IdentTypeMap.put("SPROT_NAME", "SwissProt");
        IdentTypeMap.put("DBJ", "Genbank");
        IdentTypeMap.put("SP", "SwissProtAccn");
        IdentTypeMap.put("IPI", "IPI");
        IdentTypeMap.put("COG", "COG");
        IdentTypeMap.put("ENSEMBL", "ENSEMBL");
        IdentTypeMap.put("REFSEQ_NP", "REFSEQ");
        IdentTypeMap.put("PDB", "PDB");
        IdentTypeMap.put("UNIPROT/TREMBL", "SwissProtAccn");
        IdentTypeMap.put("TREMBL", "SwissProtAccn");
        IdentTypeMap.put("REFSEQ_XP", "REFSEQ");
        IdentTypeMap.put("ORFP", "SGD_LOCUS");
        IdentTypeMap.put("UNIPROT/SPROT", "SwissProtAccn");
        IdentTypeMap.put("SWISS-PROT", "SwissProtAccn");
        IdentTypeMap.put("TPG", "Genbank");
        IdentTypeMap.put("UG", "Unigene");
        IdentTypeMap.put("SI", "SI");
        IdentTypeMap.put("UPTR", "SwissProtAccn");
        IdentTypeMap.put("UPSP", "SwissProt");
        IdentTypeMap.put("GP", "Genbank");
        IdentTypeMap.put("PIR", "PIR");
        IdentTypeMap.put("PIR2", "PIR");
    }



    public Protein(String header, byte[] bytes)
    {
        m_bytes = bytes;
        int firstAliasIndex = 0;

        m_origHeader = header;

        // Check for special case of repeated gi| at start... if so, remove the initial text, but use it for lookup string
        if (header.startsWith("gi|"))
            {
            firstAliasIndex = header.indexOf(" gi|", 2) + 1;
            if (firstAliasIndex < 0 || firstAliasIndex > 30)
                firstAliasIndex = 0;
            }

        if (0 == firstAliasIndex)
            {
            header = header.replaceAll("\t", " "); // Some annoying FASTA files have tabs instead of spaces 

            int firstSpace = header.indexOf(" ");

            if (-1 != firstSpace)
                m_lookup = header.substring(0, firstSpace).trim();
            else
                m_lookup = header;

            if (m_lookup.length() > 79)
                m_lookup = m_lookup.substring(0, 79);   // Comet truncates uni.tubingen.protein after first 79 characters
            }
        else
            m_lookup = header.substring(0, firstAliasIndex).trim();

        int massStart = header.lastIndexOf("[MASS=");

        if (massStart >= 0)
        {
            try
                {
                int massEnd = header.indexOf(']', massStart);
                m_mass = Double.parseDouble(header.substring(massStart + 6, massEnd));
                }
            catch(Exception e)
                {
                // fall through
                }
        }
        else
            massStart = header.length();

        if (0 == m_mass)
            m_mass = PeptideGenerator.computeMass(m_bytes, 0, m_bytes.length, PeptideGenerator.AMINO_ACID_AVERAGE_MASSES);

        m_header = header.substring(firstAliasIndex, massStart);
    }

    public String getHeader()
    {
        return m_header;
    }

    public void setHeader(String header)
    {
        m_header = header;
    }

    public String getOrigHeader() {
       return m_origHeader;
    }
    public void setOrigHeader(String h) {
       this.m_origHeader = h;
    }

    public String getName()
    {
        return getHeader().substring(0, Math.min(getHeader().length(), 80));
    }

    public String toString()
    {
        return getName();
    }

    public byte[] getBytes()
    {
        return m_bytes;
    }

    public void setBytes(byte[] bytes)
    {
        m_bytes = bytes;
    }

    public String getSequenceAsString()
    {
        return new String(getBytes());
    }

    public Alias[] getAliases()
        {
        String[] aliasStrings = m_header.split("\01");
        Alias[] aliases = new Alias[aliasStrings.length];

        for (int i=0; i<aliasStrings.length; i++)
            aliases[i] = new Alias(aliasStrings[i]);

        return aliases;
        }

    public double getMass()
        {
        return m_mass;
        }

    public String getLookup()
        {
        return m_lookup;
        }

    //lazily parse the header for identifiers
    public Map getIdentifierMap()
    {
        if (m_identifierMap == null)
        {
            String lookupString = m_lookup;
            if (lookupString.startsWith("IPI") && !lookupString.contains("|") && m_header.contains(" "))
            {
                lookupString = m_header.substring(m_header.indexOf(" ") + 1);
            }
            m_identifierMap = identParse(lookupString);
        }
        return m_identifierMap;
    }

    /**
     * Save out to a PrintWriter in fasta format
     * @param out
     */
    public void saveFastaFormat(PrintWriter out)
    {
        out.println(">" + m_header);
        out.println(new String(m_bytes));
        out.flush();
    }

    /**
     * Save a uni.tubingen.protein array in fasta format
     * @param proteins
     * @param outFastaFile
     * @return
     */
    public static void saveProteinArrayToFasta(ArrayList proteins,
                                               File outFastaFile)
            throws IOException
    {
        PrintWriter pw = null;
        try
            {
            pw = new PrintWriter(new FileOutputStream(outFastaFile));
            saveProteinArrayToFasta(proteins,pw);
            }
        catch (IOException x)
            {
                throw x;
            }
           finally
                {
                if (null != pw)
                    pw.close();
                }
    }

    /**
     * Save a uni.tubingen.protein array in fasta format
     * @param proteins
     * @param pw
     */
    public static void saveProteinArrayToFasta(ArrayList proteins, PrintWriter pw)
    {
        for (int i=0; i<proteins.size(); i++)
        {
            ((Protein)proteins.get(i)).saveFastaFormat(pw);
        }
    }

    /**
     *
     * @param s
     * @return
     */
    public static boolean mightBeASwissProtName(String s)
    {
        if (s.indexOf("_") == -1) return false;
        if (s.indexOf(":") != -1) s = s.substring(s.indexOf(":") + 1);
        String halves[] = s.split("_");
        if (halves.length != 2) return false;
        if (halves[0].length() < 1 || halves[0].length() > 6) return false;
        return !(halves[1].length() < 3 || halves[1].length() > 5);
    }


    /**
     *  Do some extensive parsing on a Fasta Header
     *  @return a Map, the keys are ident-types.  the values are a Set of
     *  identifier values for that type
     */
    public static Map identParse(String fastaIdentifierString)
    {
        Map identifiers = new HashMap();
        if (fastaIdentifierString.indexOf(" ") != -1) fastaIdentifierString = fastaIdentifierString.substring(0, fastaIdentifierString.indexOf(" "));
        fastaIdentifierString = fastaIdentifierString.replaceAll(":", "|");
        fastaIdentifierString = fastaIdentifierString.replace("|$", "");
        fastaIdentifierString = fastaIdentifierString.replace("|COG", "|COG|COG");

        if (fastaIdentifierString.indexOf(SEPARATOR_CHAR) == -1)
        {
            if (fastaIdentifierString.startsWith("IPI"))
            {
                fastaIdentifierString = "IPI|" + fastaIdentifierString;
            }
            else
            {
                return identifiers;
            }
        }
        String tokens[] = fastaIdentifierString.split(SEPARATOR_PATTERN);
        for (int i = 0; i < tokens.length; i++)
        {
            String key = tokens[i];
            if (key.equalsIgnoreCase("gnl"))
            {
                if (i < (tokens.length - 2))
                {
                    key = tokens[++i];
                }
            }
            String value = null;
            if (i + 1 < tokens.length)
            {
                value = tokens[++i];
            }
            if (!(tokens.length == 2 & tokens[0].equalsIgnoreCase("UPSP")) && i == (tokens.length - 1) && mightBeASwissProtName(tokens[i]))
            {
                value = key;
            }
            if (key.equalsIgnoreCase("SI") && value != null)
            {
                String halves[] = value.split("_");
                value = halves[0];
            }
            // if we really found an identifier, put it into the set of identifiers
            // of the same type.  Make sure it isn't too big (50)

            if (value != null && IdentTypeMap.containsKey(key.toUpperCase()))
            {
                String newKey = (String)IdentTypeMap.get(key.toUpperCase());
                if (!identifiers.containsKey(newKey))
                    identifiers.put(newKey, new HashSet());
                Set curSet = (Set)identifiers.get(newKey);
                String commaSepGroup[] = value.split(";");
                for (int k=0;k<commaSepGroup.length;k++)
                {
                	String oneOfList  = commaSepGroup[k];
                    if (oneOfList.length() > 50) oneOfList = oneOfList.substring(0, 49);
                    curSet.add(oneOfList);
                }
            }
        }
        return identifiers;
    }





    public class Alias
        {
        private String _remaining = "";
        private String _description = "";

        public Alias(String s)
            {
            _description = s;
            }

        public String getRemaining()
            {
            return _remaining;
            }

        public String getDescription()
            {
            return _description;
            }
        }

    public static class SequenceComparator implements Comparator
        {
        public int compare(Object o1, Object o2)
            {
            return ((Protein) o1).getSequenceAsString().compareTo(((Protein) o2).getSequenceAsString());

            }
        }

    public static class LookupComparator implements Comparator
        {
        public int compare(Object o1, Object o2)
            {
            return ((Protein) o1).getLookup().compareTo(((Protein) o2).getLookup());

            }
        }

    public static class HeaderComparator implements Comparator
        {
        public int compare(Object o1, Object o2)
            {
            return ((Protein) o1).getHeader().compareTo(((Protein) o2).getHeader());

            }
        }
}
