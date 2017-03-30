package bioproject;

import static bioproject.MainPage.three;
import java.awt.Color;
import java.awt.Toolkit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.JOptionPane;
import javax.swing.text.BadLocationException;
import javax.swing.text.DefaultHighlighter;
import javax.swing.text.Highlighter;

/**
 *
 * @author Nize
 */
public class Three extends javax.swing.JFrame {

    int count = 0;
    Algorithms a = new Algorithms();
    String txt;

    public Three() {
        initComponents();
        getIcon();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jLabel2 = new javax.swing.JLabel();
        jLabel3 = new javax.swing.JLabel();
        jLabel4 = new javax.swing.JLabel();
        jScrollPane1 = new javax.swing.JScrollPane();
        jTextArea1 = new javax.swing.JTextArea();
        jScrollPane2 = new javax.swing.JScrollPane();
        jTextArea2 = new javax.swing.JTextArea();
        jLabel5 = new javax.swing.JLabel();
        jLabel6 = new javax.swing.JLabel();
        jLabel7 = new javax.swing.JLabel();
        jLabel1 = new javax.swing.JLabel();

        setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);
        setTitle("หา Gene prediction");
        setMaximumSize(new java.awt.Dimension(732, 499));
        setMinimumSize(new java.awt.Dimension(732, 499));
        setPreferredSize(new java.awt.Dimension(732, 530));
        setResizable(false);
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                formWindowClosing(evt);
            }
        });
        getContentPane().setLayout(null);

        jLabel2.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jLabel2MouseClicked(evt);
            }
        });
        getContentPane().add(jLabel2);
        jLabel2.setBounds(690, 160, 40, 40);

        jLabel3.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jLabel3MouseClicked(evt);
            }
        });
        getContentPane().add(jLabel3);
        jLabel3.setBounds(490, 230, 190, 40);

        jLabel4.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jLabel4MouseClicked(evt);
            }
        });
        getContentPane().add(jLabel4);
        jLabel4.setBounds(310, 460, 130, 30);

        jTextArea1.setColumns(20);
        jTextArea1.setFont(new java.awt.Font("Monospaced", 0, 16)); // NOI18N
        jTextArea1.setRows(5);
        jTextArea1.setText("กรอกสาย DNA หรือ mRNA ที่ต้องการแปลง");
        jTextArea1.addCaretListener(new javax.swing.event.CaretListener() {
            public void caretUpdate(javax.swing.event.CaretEvent evt) {
                jTextArea1CaretUpdate(evt);
            }
        });
        jTextArea1.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jTextArea1MouseClicked(evt);
            }
        });
        jTextArea1.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent evt) {
                jTextArea1KeyTyped(evt);
            }
        });
        jScrollPane1.setViewportView(jTextArea1);

        getContentPane().add(jScrollPane1);
        jScrollPane1.setBounds(50, 40, 640, 160);

        jTextArea2.setEditable(false);
        jTextArea2.setColumns(20);
        jTextArea2.setFont(new java.awt.Font("Monospaced", 0, 16)); // NOI18N
        jTextArea2.setRows(5);
        jScrollPane2.setViewportView(jTextArea2);

        getContentPane().add(jScrollPane2);
        jScrollPane2.setBounds(50, 290, 640, 160);

        jLabel5.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jLabel5MouseClicked(evt);
            }
        });
        getContentPane().add(jLabel5);
        jLabel5.setBounds(70, 10, 110, 30);

        jLabel6.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jLabel6MouseClicked(evt);
            }
        });
        getContentPane().add(jLabel6);
        jLabel6.setBounds(690, 10, 40, 40);

        jLabel7.setFont(new java.awt.Font("Tahoma", 0, 12)); // NOI18N
        jLabel7.setForeground(new java.awt.Color(153, 0, 0));
        jLabel7.setText("หาด้วยตัวเอง");
        jLabel7.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jLabel7MouseClicked(evt);
            }
        });
        getContentPane().add(jLabel7);
        jLabel7.setBounds(630, 460, 70, 20);

        jLabel1.setIcon(new javax.swing.ImageIcon(getClass().getResource("/Photo/03.jpg"))); // NOI18N
        getContentPane().add(jLabel1);
        jLabel1.setBounds(0, 0, 940, 500);

        pack();
        setLocationRelativeTo(null);
    }// </editor-fold>//GEN-END:initComponents

    private void formWindowClosing(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosing
        MainPage m = new MainPage();
        m.setVisible(true);
        this.setVisible(false);
    }//GEN-LAST:event_formWindowClosing

    private void jLabel2MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jLabel2MouseClicked
        jTextArea1.setText("กรอกสาย DNA หรือ mRNA ที่ต้องการแปลง");
        count = 0;
    }//GEN-LAST:event_jLabel2MouseClicked

    private void jTextArea1KeyTyped(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_jTextArea1KeyTyped
        if (count == 0) {
            jTextArea1.setText("");
        }
        count++;
    }//GEN-LAST:event_jTextArea1KeyTyped

    private void jTextArea1MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jTextArea1MouseClicked
        if (count == 0) {
            jTextArea1.setText("");
        }
        count++;
    }//GEN-LAST:event_jTextArea1MouseClicked

    private void jLabel3MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jLabel3MouseClicked
        txt = jTextArea1.getText();
        if (txt.length() > 0 && txt.charAt(0) != 'ก') {
            if (!a.CheckError(txt.replaceAll("\n", ""))) {
                if (a.GenePrediction(a.DNATomRNA(txt.toUpperCase().replaceAll("\n", ""))).length() > 0) {
                    jTextArea2.setText(a.GenePrediction(a.DNATomRNA(txt.toUpperCase().replaceAll("\n", ""))));
                } else {
                    jTextArea2.setText("N/A");
                    JOptionPane.showMessageDialog(this,
                            "ไม่พบรหัส codon ใดที่สามารถแปลงเป็นกรดอมิโน", "แจ้งเตือน",
                            JOptionPane.ERROR_MESSAGE);
                }

            } else if (JOptionPane.showConfirmDialog(null, "สายข้อมูลนำเข้ามีข้อผิดพลาด\nต้องการดำเนินการต่อหรือไม่", "แจ้งเตือน",
                    JOptionPane.YES_NO_OPTION) == JOptionPane.YES_OPTION) {
                if (a.GenePrediction(a.DNATomRNA(txt.toUpperCase().replaceAll("\n", ""))).length() > 0) {
                    jTextArea2.setText(a.GenePrediction(a.DNATomRNA(txt.toUpperCase().replaceAll("\n", ""))));
                } else {
                    jTextArea2.setText("N/A");
                    JOptionPane.showMessageDialog(this,
                            "ไม่พบรหัส codon ใดที่สามารถแปลงเป็นกรดอมิโน", "แจ้งเตือน",
                            JOptionPane.ERROR_MESSAGE);
                }
            } else {
                jTextArea2.setText("N/A");
            }
        } else {
            JOptionPane.showMessageDialog(this,
                    "กรุณากรอกสาย DNA หรือ mRNA ที่ต้องการแปลง", "แจ้งเตือน",
                    JOptionPane.WARNING_MESSAGE);
        }
    }//GEN-LAST:event_jLabel3MouseClicked
    public void getCo(String a) {
        Highlighter.HighlightPainter cyanPainter = new DefaultHighlighter.DefaultHighlightPainter(Color.cyan);
        Highlighter.HighlightPainter redPainter = new DefaultHighlighter.DefaultHighlightPainter(Color.red);
        String tmp = a;
        Pattern p = Pattern.compile("ATG");
        Matcher m = p.matcher(tmp);
        while (m.find()) {
            try {
                jTextArea1.getHighlighter().addHighlight(m.start(), m.end(), cyanPainter);
            } catch (BadLocationException ex) {
                System.out.println("Not found.");
            }
        }
        p = Pattern.compile("TTA");
        m = p.matcher(tmp);
        while (m.find()) {
            try {
                jTextArea1.getHighlighter().addHighlight(m.start(), m.end(), redPainter);
            } catch (BadLocationException ex) {
                System.out.println("Not found.");
            }
        }
        p = Pattern.compile("TAG");
        m = p.matcher(tmp);
        while (m.find()) {
            try {
                jTextArea1.getHighlighter().addHighlight(m.start(), m.end(), redPainter);
            } catch (BadLocationException ex) {
                System.out.println("Not found.");
            }
        }
        p = Pattern.compile("TGA");
        m = p.matcher(tmp);
        while (m.find()) {
            try {
                jTextArea1.getHighlighter().addHighlight(m.start(), m.end(), redPainter);
            } catch (BadLocationException ex) {
                System.out.println("Not found.");
            }
        }
    }
    public void getCo1(String a) {
        Highlighter.HighlightPainter cyanPainter = new DefaultHighlighter.DefaultHighlightPainter(Color.cyan);
        Highlighter.HighlightPainter redPainter = new DefaultHighlighter.DefaultHighlightPainter(Color.red);
        Pattern p = Pattern.compile("[^atgcuATGCU\n]");
        Matcher m = p.matcher(a);
        while (m.find()) {
            try {
                jTextArea1.getHighlighter().addHighlight(m.start(), m.end(), redPainter);
            } catch (BadLocationException ex) {
                System.out.println("Not found.");
            }
        }
    }
    private void jLabel4MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jLabel4MouseClicked
        jTextArea2.setText("");
    }//GEN-LAST:event_jLabel4MouseClicked

    private void jLabel5MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jLabel5MouseClicked
        jTextArea1.setText("ATGCGCTCTCTTTCAGCATTCGTCTTCCTCTCGGTCATCGTCGGCACA\n"
                + "TATGGCGCAATAGGTCCAAGCACCAGTCTTTACATCGCAAACAAGTAC\n"
                + "ATCCTCCCCGATGGATTCAACCGCTCTAGTGTTTTGGCAGGTCCCACT\n"
                + "GCAGCAAATGTATCATTCCCTGGTCCTGTTATCACTGGCTTCCAGGGA\n"
                + "GATACATTTAGCATTAATGTCATCGATCAACTCACTGATACTACAATG\n"
                + "TTGACGAGCACTAGTCTTCACTGGCATGGCCTGTTCTAAGAAGGTAGT");
    }//GEN-LAST:event_jLabel5MouseClicked

    private void jLabel6MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jLabel6MouseClicked
        About03 ab = new About03();
        ab.setVisible(true);
        three.setEnabled(false);
    }//GEN-LAST:event_jLabel6MouseClicked

    private void jLabel7MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jLabel7MouseClicked
        if (jTextArea1.getText().length() > 0 && jTextArea1.getText().charAt(0) != 'ก') {
            getCo(jTextArea1.getText());
        } else {
            JOptionPane.showMessageDialog(this,
                    "กรุณากรอกสาย DNA หรือ mRNA ที่ต้องการแปลง", "แจ้งเตือน",
                    JOptionPane.WARNING_MESSAGE);
        }

    }//GEN-LAST:event_jLabel7MouseClicked

    private void jTextArea1CaretUpdate(javax.swing.event.CaretEvent evt) {//GEN-FIRST:event_jTextArea1CaretUpdate
        getCo1(jTextArea1.getText());
    }//GEN-LAST:event_jTextArea1CaretUpdate

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(Three.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(() -> {
            new Three().setVisible(true);
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JTextArea jTextArea1;
    private javax.swing.JTextArea jTextArea2;
    // End of variables declaration//GEN-END:variables

    private void getIcon() {
        setIconImage(Toolkit.getDefaultToolkit().getImage(getClass().getResource("icon.ico")));
    }
}
