package uni.tubingen.inference.lp.algorithm.protein;

import java.io.*;
import java.util.*;

import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.FileInputStream;
import java.io.InputStream;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.ProgressMonitorInputStream;

/**
 * Date: Sep 14, 2007
 * Time: 20:17:58 PM
 *
 */
public class FastaLoader
{
    private File m_fastaFile;

    public FastaLoader(File fastaFile)
    {
        m_fastaFile = fastaFile;
    }

    public Iterator iterator()
    {
        return new ProteinIterator();
    }

    public class ProteinIterator implements Iterator
    {
        String m_proteinHeader = null;
        BufferedReader m_reader = null;
        private boolean m_beforeFirst = true;

        private void init()
        {
            try
            {
                m_reader = new BufferedReader(new FileReader(m_fastaFile));
                String line = m_reader.readLine();

                //Iterator expects m_proteinHeader to be initialized...
                if (null != line && line.charAt(0) == '>')
                    m_proteinHeader = line.substring(1);
                else
                {
                    if (null != m_reader)
                        m_reader.close();

                    throw new IllegalArgumentException("Fasta File did not start with a >");
                }
            }
            catch (IOException x)
            {
                if (null != m_reader)
                    try
                    {
                        m_reader.close();
                    }
                    catch (IOException x2) {}

                throw new IllegalArgumentException ("Error reading from " + m_fastaFile.getAbsolutePath());
            }

            m_beforeFirst = false;
        }

        /**
         *
         * @return are there any more proteins left in the file
         */
        public boolean hasNext()
        {
            if (m_beforeFirst)
                init();

            return null != m_proteinHeader;
        }

        /**
         * Closes file just in case.
         * @throws IOException if file is not closeable
         */
        protected void finalize() throws Throwable
        {
            super.finalize();    //If iteration is not complete, still close the file...
            if (null != m_reader)
                m_reader.close();
        }

        /**
         * Get next uni.tubingen.protein object in file.
         * @return Protein or null if end of file
         */
        public Object next()
        {
            if (m_beforeFirst)
                init();

            if (null == m_proteinHeader)
                return null;

            ByteArrayOutputStream aaStream = new ByteArrayOutputStream(2048);
            String line;

            try
            {
                while((line = m_reader.readLine()) != null)
                {
                    if (line.length() > 0 && line.charAt(0) == '>')
                    {
                        Protein p = new Protein(m_proteinHeader, aaStream.toByteArray());
                        m_proteinHeader = line.substring(1);
                        return p;
                    }
                    else
                    {
                        byte[] bytes = line.getBytes();
                        for (int i = 0; i < bytes.length; i++)
                            if((bytes[i] >= 'A') && (bytes[i] <= 'Z'))
                            {
                                //_aaCounts[bytes[i] - 'A'] ++;
                                aaStream.write(bytes[i]);
                            }
                    }
                }
                
                //End of file
                Protein p = new Protein(m_proteinHeader, aaStream.toByteArray());
                close();
                return p;
            }
            catch (IOException x)
            {
                return null;
            }

        }

        /**
         * Unsupported
         */
        public void remove()
        {
            throw new UnsupportedOperationException();
        }

        /**
         * Closes the file. No more items will be returned from the iterator
         */
        public void close()
        {
            if (null != m_reader)
                try
                {
                    m_reader.close();
                }
                catch (IOException x) {}

            m_reader = null;
            m_proteinHeader = null;
        }
    }
    
    public static void main(String[] args)
    {
    	
        FastaLoader loader = null;
        try
        {
            File fastaFile = new File("Sigma_49.fasta");
            loader = new FastaLoader(fastaFile);
        }
        catch (Exception x)
        {
            x.printStackTrace(System.err);
            return;
        }

        Iterator it = loader.iterator();

        int uu=0,i=0;
        while (it.hasNext())
        {
            Protein p = (Protein) it.next();
            byte[] bytes = p.getBytes();
            System.out.println(p.getHeader());
            if(p.getHeader().startsWith("[Contaminant]")){
            	i++;
            }
            uu++;
        }
        System.out.println(uu);
        System.out.println(i);
  }
    
}
