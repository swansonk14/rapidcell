import java.awt.GridBagLayout;
import javax.swing.JPanel;
import javax.swing.JButton;
import java.awt.Rectangle;
import java.awt.event.KeyEvent;

import javax.swing.*;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.Point;
import java.awt.Toolkit;
import java.io.File;
import java.awt.FlowLayout;
import java.awt.Color;
import javax.swing.JTextArea;
import java.awt.Canvas;
import java.awt.event.KeyAdapter;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.TextArea;
import java.awt.Component;
import javax.swing.border.BevelBorder;
import java.awt.Font;

public class Panel1 extends JPanel {

	private static final long serialVersionUID = 1L;
	private JButton jButtonRun = null;
	private JButton jButtonStop = null;
	public Thread ThModel1 = null;
	private JButton jButtonOK = null;
	private JTabbedPane jTabbedPaneNetwork = null;
	private JPanel jPanelMain = null;
	private JPanel jPanelGradient = null;
	private JPanel jPanelAdvanced = null;
	private JLabel jLabel1 = null;
	private JTextField jTextFieldNcells = null;
	private JLabel jLabel11 = null;
	private JTextField jTextFieldTmax = null;
	private JLabel jLabel12 = null;
	private JTextField jTextFieldLx = null;
	private JLabel jLabel13 = null;
	private JComboBox jComboBoxIniPosition = null;
	private JLabel jLabel14 = null;
	private JComboBox jComboBoxGradient = null;
	private JComboBox jComboBoxLigandMethylation = null;
	private JLabel jLabel15 = null;
	private JComboBox jComboBoxMedium = null;
	private JLabel jLabel16 = null;
	private JTextField jTextFieldK_offTar = null;
	private JTextField jTextFieldK_onTar = null;
	private JTextField jTextFieldK_offTsr = null;
	private JTextField jTextFieldK_onTsr = null;
	private JLabel jLabel17 = null;
	private JLabel jLabel18 = null;
	private JLabel jLabel19 = null;
	private JLabel jLabel110 = null;
	private JLabel jLabel1101 = null;
	private JLabel jLabel1102 = null;
	private JLabel jLabel11011 = null;
	private JLabel jLabel111 = null;
	private JComboBox jComboBoxModelRun = null;
	private JLabel jLabel1111 = null;
	private JComboBox jComboBoxModelTumble = null;
	private JLabel jLabel11111 = null;
	private JTextField jTextFieldAdRate = null;
	private JLabel jLabel111111 = null;
	private JTextField jTextFieldNmotors = null;
	private JLabel jLabel111112 = null;
	private JTextField jTextFieldDr = null;
	private JLabel jLabel1111111 = null;
	private JTextField jTextFieldVswim = null;
	private JLabel jLabel1111112 = null;
	private JTextField jTextFieldNtar = null;
	private JTextField jTextFieldNtsr = null;
	private JLabel jLabel112 = null;
	private JLabel jLabel1121 = null;
	private JLabel jLabel11111121 = null;
	private JTextField jTextFieldMbias0 = null;
	private JProgressBar jProgressBar = null;
	private JPanel jPanelNetwork = null;
	private JLabel jLabel113 = null;
	private JTextField jTextFieldAdPrecision = null;
	private JLabel jLabel114 = null;
	private JLabel jLabel1141 = null;
	private JTextField jTextFieldCheR = null;
	private JTextField jTextFieldCheB = null;
	private JComboBox jComboBoxModelMotor = null;
	private JLabel jLabel115 = null;
	private JTextField jTextFieldCheY = null;
	private JTextField jTextFieldCheZ = null;
	private JLabel jLabel116 = null;
	private JLabel jLabel117 = null;
	private JTextField jTextFieldMeth = null;
	private JLabel jLabel120 = null;
	private JTextField jTextFieldStep = null;
	private JLabel jLabel121 = null;
	private JTextField jTextFieldAdditionTime = null;
	private JLabel jLabel122 = null;
	private JTextField jTextFieldRemovalTime = null;
	private JLabel jLabel123 = null;
	private JLabel jLabelGauss1 = null;
	private JTextField jTextFieldMaxAsp = null;
	private JLabel jLabelGradMin = null;
	private JTextField jTextFieldMinAsp = null;
	private JLabel jLabelGradCenterXY = null;
	private JTextField jTextFieldLy = null;
	private JTextField jTextFieldGradCenterX = null;
	private JTextField jTextFieldGradCenterY = null;
	private JLabel jLabelGradRate = null;
	private JTextField jTextFieldGradRate = null;
	private JLabel jLabelSigma = null;
	private JTextField jTextFieldSigma = null;
	private JPanel jPanelGradParam2 = null;
	private JPanel jPanelGradParam1 = null;
	private JPanel jPanelKd = null;
	private JTabbedPane jTabbedPaneOutput = null;
	private JTextField jTextFieldIniPosX = null;
	private JTextField jTextFieldIniPosY = null;
	private JLabel jLabelIniOri = null;
	private JComboBox jComboBoxIniOri = null;
	private JPanel jPanelOutput = null;
	private JLabel jLabelOutDt = null;
	private JTextField jTextFieldOutDt = null;
	private JLabel jLabelOutFileName1 = null;
	private JLabel jLabelOutFileName2 = null;
	private JTextField jTextFieldFileNameAve = null;
	private JTextField jTextFieldFileNameFull = null;
	private JLabel jLabelOutXY = null;
	private JCheckBox jCheckBoxOutXY = null;
	private JLabel jLabelOutOri = null;
	private JCheckBox jCheckBoxOutOri = null;
	private JLabel jLabelOutCheA = null;
	private JCheckBox jCheckBoxOutCheA = null;
	private JLabel jLabelOutCheY = null;
	private JCheckBox jCheckBoxOutCheY = null;
	private JLabel jLabelOutMeth = null;
	private JCheckBox jCheckBoxOutMeth = null;
	private JLabel jLabelOutMb = null;
	private JCheckBox jCheckBoxOutMb = null;
	private JLabel jLabelOutRun = null;
	private JCheckBox jCheckBoxOutRun = null;
	private JLabel jLabelOutAttr = null;
	private JCheckBox jCheckBoxOutAttr = null;
	private JLabel jLabelOutNmot = null;
	private JCheckBox jCheckBoxOutNmot = null;
	private JLabel jLabelOutSwitchFreq = null;
	private JCheckBox jCheckBoxOutSwitchFreq = null;
	private JCheckBox jCheckBoxOutPon = null;
	private JLabel jLabelCheckAll = null;
	private JCheckBox jCheckBoxCheckAll = null;
	private JLabel jLabelUncheckAll = null;
	private JCheckBox jCheckBoxUncheckAll = null;
	private JPanel jPanelAlgorithm = null;
	private JLabel jLabelDt = null;
	private JTextField jTextFieldDt = null;
	private JLabel jLabelSeed = null;
	private JTextField jTextFieldSeed = null;
	private JLabel jLabelBoundary = null;
	private JComboBox jComboBoxBoundary = null;
	private ArenaCanvas canvas;
	private JTextArea textAreaConsole;
	private JPanel panelArena;
	private JTextField jTextFieldCheA;
	private JPanel jpanelCluster;
	private JTextField textFieldkA;
	private JTextField textFieldkY;
	private JTextField textFieldkZ;
	private JTextField textFieldkR;
	private JTextField textFieldkB;
	private JPanel panelCheckBoxes;
	/**
	 * This is the default constructor
	 */
	public Panel1() {
		super();
		setSize(new Dimension(785, 396));
		initialize();
	}

	/**
	 * This method initializes this
	 * 
	 * @return void
	 */
	private void initialize() {
		this.setLayout(null);
		this.add(getJButtonRun(), null);
		this.add(getJButtonStop(), null);
		this.add(getJButtonOK(), null);
		this.add(getJTabbedPaneNetwork(), null);
		this.add(getJProgressBar(), null);
		
		textAreaConsole = new JTextArea("Welcome to RapidCell 1.4");
		this.add(textAreaConsole);
		textAreaConsole.setFont(new Font("Tahoma", Font.PLAIN, 14));
		textAreaConsole.setEditable(false);

		JScrollPane scrollPane = new JScrollPane(textAreaConsole);
		scrollPane.setBounds(12, 400, 760, 130);
		add(scrollPane);

		setJComboBoxExponentialDisplacementGradient();
	}
	
	public void resetPanel(){
		jButtonRun.setVisible(true);
		jButtonStop.setVisible(false);
	}
	public void displayResults(){
		File directory=new File (".");
		jButtonStop.setVisible(false);
		textAreaConsole.append("\n Done! Output is written to folder: \n");
		textAreaConsole.setCaretPosition( textAreaConsole.getDocument().getLength() );
		try {textAreaConsole.append(directory.getCanonicalPath());}
			catch(Exception e) {textAreaConsole.append("\n Error: Exception is ="+e.getMessage());}
		jButtonOK.setVisible(true);
	}
	public void hideResults(){
		jButtonOK.setVisible(false);
		resetPanel();
		canvas.restart();
	}
	public void readParams(Model myModel){
		myModel.setNcells(Integer.parseInt(jTextFieldNcells.getText()));
		myModel.setTmax(Double.parseDouble(jTextFieldTmax.getText()));
		myModel.setLxLy(Double.parseDouble(jTextFieldLx.getText()),Double.parseDouble(jTextFieldLy.getText()));
		myModel.setInitialCellDistribution(Double.parseDouble(jTextFieldIniPosX.getText()),
				Double.parseDouble(jTextFieldIniPosY.getText()),jComboBoxIniOri.getSelectedIndex());
		myModel.setGradient(jComboBoxGradient.getSelectedIndex());
		myModel.setStepGradientParams(Double.parseDouble(jTextFieldStep.getText()),
				Double.parseDouble(jTextFieldAdditionTime.getText()),
				Double.parseDouble(jTextFieldRemovalTime.getText()));
		myModel.setGradMinMaxVal(Double.parseDouble(jTextFieldMinAsp.getText()), Double.parseDouble(jTextFieldMaxAsp.getText()));
		myModel.setGradSourceXY(Double.parseDouble(jTextFieldGradCenterX.getText()), Double.parseDouble(jTextFieldGradCenterY.getText()));
		myModel.setExpGradRate(Double.parseDouble(jTextFieldGradRate.getText()));
		myModel.setGradSigma(Double.parseDouble(jTextFieldSigma.getText()));
		myModel.setMedium(jComboBoxMedium.getSelectedIndex());
		myModel.setOutputDt(Double.parseDouble(jTextFieldOutDt.getText()));
		myModel.setOutFileNames(jTextFieldFileNameAve.getText(), jTextFieldFileNameFull.getText());
		myModel.setOutFlags(jCheckBoxOutXY.isSelected(), jCheckBoxOutXY.isSelected(), 
				jCheckBoxOutOri.isSelected(), jCheckBoxOutCheY.isSelected(), jCheckBoxOutCheA.isSelected(), 
				jCheckBoxOutMeth.isSelected(), jCheckBoxOutMb.isSelected(), 
				jCheckBoxOutSwitchFreq.isSelected(), jCheckBoxOutNmot.isSelected(), 
				jCheckBoxOutRun.isSelected(), jCheckBoxOutAttr.isSelected(), jCheckBoxOutPon.isSelected());
		myModel.setDt(Double.parseDouble(jTextFieldDt.getText()));
		myModel.setRandSeed(Long.parseLong(jTextFieldSeed.getText()));
		Cell2D.setRandSeed(Long.parseLong(jTextFieldSeed.getText()));
		Network.setRandSeed(Long.parseLong(jTextFieldSeed.getText()));
		Network.setDissConstants(Double.parseDouble(jTextFieldK_offTar.getText()), 
				Double.parseDouble(jTextFieldK_onTar.getText()), 
				Double.parseDouble(jTextFieldK_offTsr.getText()),
				Double.parseDouble(jTextFieldK_onTsr.getText()));
		myModel.setRunModel(jComboBoxModelRun.getSelectedIndex());
		myModel.setTumbleModel(jComboBoxModelTumble.getSelectedIndex());
		myModel.setBoundaryCondition(jComboBoxBoundary.getSelectedIndex());
		Cell.setNmotors(Integer.parseInt(jTextFieldNmotors.getText()));
		Cell.setRotDiffusionCoeff(Double.parseDouble(jTextFieldDr.getText()));
		Cell.setMaxCellVelocity(1e-3*Double.parseDouble(jTextFieldVswim.getText()));
		Motor.setMotorBias0(Double.parseDouble(jTextFieldMbias0.getText()),jComboBoxModelMotor.getSelectedIndex());
		//Panel Network:
		myModel.setAdaptRate(Double.parseDouble(jTextFieldAdRate.getText()));
		Network.setAdaptPrecision(0.01*Double.parseDouble(jTextFieldAdPrecision.getText()));
		Network.setRelativeCheRCheB(0.01*Double.parseDouble(jTextFieldCheR.getText()), 0.01*Double.parseDouble(jTextFieldCheB.getText()));
		Network.setNreceptors(Integer.parseInt(jTextFieldNtar.getText()), Integer.parseInt(jTextFieldNtsr.getText()));
		Network.setRelativeCheYCheZ(0.01*Double.parseDouble(jTextFieldCheY.getText()), 0.01*Double.parseDouble(jTextFieldCheZ.getText()));
		Network.setIniMethState(Double.parseDouble(jTextFieldMeth.getText()));
		Network.setRelativeCheA(0.01*Double.parseDouble(jTextFieldCheA.getText()));
		Network.setRates(Double.parseDouble(textFieldkR.getText()), Double.parseDouble(textFieldkB.getText()), 
				Double.parseDouble(textFieldkA.getText()), Double.parseDouble(textFieldkY.getText()), Double.parseDouble(textFieldkZ.getText()));
		Network.setLigandMethylationRelation(jComboBoxLigandMethylation.getSelectedIndex());
	}
	/**
	 * This method initializes jButtonRun	
	 * 	
	 * @return javax.swing.JButton	
	 */
	private JButton getJButtonRun() {
		if (jButtonRun == null) {
			jButtonRun = new JButton();
			jButtonRun.setFont(new Font("Tahoma", Font.PLAIN, 14));
			jButtonRun.setBounds(new Rectangle(18, 358, 113, 35));
			jButtonRun.setMnemonic(KeyEvent.VK_UNDEFINED);
			jButtonRun.setText("Run");
			jButtonRun.setName("");
			jButtonRun.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					ThModel1 = new Thread() {
						Model Model1=new Model(this);
						public void run(){
							try {
								readParams(Model1);
								textAreaConsole.append("\n running...");
								Model1.runModel();
								textAreaConsole.append("\n execution time: "+ Model1.exectime/1000 +" s");
								canvas.start();
								setCanvas();
								canvas.setCells(Model1.arrCells);
								canvas.repaint();
						    	} catch (InterruptedException e) {
						    		textAreaConsole.append("\n aborted :(");
							        return;
							    }
							    if(!this.isInterrupted()){
								    SwingUtilities.invokeLater(new Runnable() {
							              public void run() {
							            	  jProgressBar.setIndeterminate(false);
							            	  Toolkit.getDefaultToolkit().beep();
							            	  displayResults();
							              }
							            });							    	
							    }
							    else { 
							    	SwingUtilities.invokeLater(new Runnable() {
							              public void run() {
							            	  textAreaConsole.append("\n interrupted :(");
							            	  textAreaConsole.setCaretPosition( textAreaConsole.getDocument().getLength() );
							              }
							              });
							    	}

						}
					};
					ThModel1.start();
					jProgressBar.setIndeterminate(true);
					jButtonRun.setVisible(false);
					jButtonStop.setVisible(true);
				}
			});
		}
		return jButtonRun;
	}

	/**
	 * This method initializes jButtonStop	
	 * 	
	 * @return javax.swing.JButton	
	 */
	private JButton getJButtonStop() {
		if (jButtonStop == null) {
			jButtonStop = new JButton();
			jButtonStop.setFont(new Font("Tahoma", Font.PLAIN, 14));
			jButtonStop.setBounds(new Rectangle(143, 358, 106, 35));
			jButtonStop.setText("Stop");
			jButtonStop.setVisible(false);
			jButtonStop.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					ThModel1.interrupt();
					ThModel1=null;
					jProgressBar.setIndeterminate(false);
					jButtonStop.setVisible(false);
					jButtonRun.setVisible(true);
				}
			});
		}
		return jButtonStop;
	}

	/**
	 * This method initializes jButtonOK	
	 * 	
	 * @return javax.swing.JButton	
	 */
	private JButton getJButtonOK() {
		if (jButtonOK == null) {
			jButtonOK = new JButton();
			jButtonOK.setFont(new Font("Tahoma", Font.PLAIN, 14));
			jButtonOK.setBounds(new Rectangle(261, 358, 106, 35));
			jButtonOK.setText("Restart");
			jButtonOK.setVisible(false);
			jButtonOK.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					hideResults(); 
				}
			});
		}
		return jButtonOK;
	}

	/**
	 * This method initializes jTabbedPaneNetwork	
	 * 	
	 * @return javax.swing.JTabbedPane	
	 */
	private JTabbedPane getJTabbedPaneNetwork() {
		if (jTabbedPaneNetwork == null) {
			jTabbedPaneNetwork = new JTabbedPane();
			jTabbedPaneNetwork.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTabbedPaneNetwork.setBounds(new Rectangle(0, 0, 780, 350));
			jTabbedPaneNetwork.addTab("Main", null, getJPanelMain(), null);
			jTabbedPaneNetwork.addTab("Gradient", null, getJPanelGradient(), null);
			jTabbedPaneNetwork.addTab("Motility", null, getJPanelAdvanced(), null);
			jTabbedPaneNetwork.addTab("Network", null, getJPanelNetwork(), null);
			jTabbedPaneNetwork.addTab("Algorithm", null, getJPanelAlgorithm(), null);
			jPanelMain.add(getCanvas());
			jTabbedPaneNetwork.setEnabledAt(1, true);
			jTabbedPaneNetwork.addTab("Output", null, getJPanelOutput(), null);
			canvas.addMouseListener(new MouseAdapter() {
				@Override
				public void mouseEntered(MouseEvent arg0) {
					setCanvas();				
					canvas.repaint();
				}
			});	
		}
		return jTabbedPaneNetwork;
	}

	
	/**
	 * This method initializes jPanelMain	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getJPanelMain() {
		if (jPanelMain == null) {
			jLabelSeed = new JLabel();
			jLabelSeed.setFont(new Font("Tahoma", Font.PLAIN, 14));
			jLabelSeed.setBounds(new Rectangle(547, 82, 98, 20));
			jLabelSeed.setText("Random seed");
			jLabel11011 = new JLabel();
			jLabel11011.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel11011.setText("Kon");
			jLabel11011.setBounds(new Rectangle(131, 70, 31, 16));
			jLabel1102 = new JLabel();
			jLabel1102.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel1102.setText("Koff");
			jLabel1102.setBounds(new Rectangle(133, 40, 31, 16));
			jLabel1101 = new JLabel();
			jLabel1101.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel1101.setText("Kon");
			jLabel1101.setBounds(new Rectangle(28, 70, 31, 16));
			jLabel110 = new JLabel();
			jLabel110.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel110.setText("Koff");
			jLabel110.setBounds(new Rectangle(28, 40, 31, 16));
			jLabel19 = new JLabel();
			jLabel19.setFont(new Font("Tahoma", Font.BOLD, 15));
			jLabel19.setText("Tsr");
			jLabel19.setBounds(new Rectangle(163, 25, 46, 16));
			jLabel18 = new JLabel();
			jLabel18.setFont(new Font("Tahoma", Font.BOLD, 15));
			jLabel18.setText("Tar");
			jLabel18.setBounds(new Rectangle(58, 25, 46, 16));
			jLabel17 = new JLabel();
			jLabel17.setFont(new Font("Tahoma", Font.PLAIN, 16));
			jLabel17.setText("Receptor-ligand Kd-s (mM)");
			jLabel17.setBounds(new Rectangle(28, 0, 221, 20));
			jLabel12 = new JLabel();
			jLabel12.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel12.setBounds(new Rectangle(279, 83, 122, 20));
			jLabel12.setHorizontalAlignment(SwingConstants.LEADING);
			jLabel12.setText("Simulation time, s");
			jLabel11 = new JLabel();
			jLabel11.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel11.setHorizontalAlignment(SwingConstants.LEFT);
			jLabel11.setBounds(new Rectangle(279, 48, 119, 20));
			jLabel11.setText("Number of cells");
			jPanelMain = new JPanel();
			jPanelMain.setFont(new Font("Tahoma", Font.PLAIN, 14));
			jPanelMain.setLayout(null);
			jPanelMain.add(getJTextFieldNcells(), null);
			jPanelMain.add(jLabel11, null);
			jPanelMain.add(getJTextFieldTmax(), null);
			jPanelMain.add(jLabel12, null);
			jPanelMain.add(jLabelSeed, null);
			jPanelMain.add(getJTextFieldSeed(), null);
			jPanelMain.add(getPanelArena());
			JPanel panelOutFiles = new JPanel();
			panelOutFiles.setBounds(279, 128, 462, 65);
			jPanelMain.add(panelOutFiles);
			panelOutFiles.setLayout(null);
			jLabelOutFileName1 = new JLabel();
			jLabelOutFileName1.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelOutFileName1.setBounds(0, 0, 163, 25);
			panelOutFiles.add(jLabelOutFileName1);
			jLabelOutFileName1.setText("Output file name (ave.)");
			jLabelOutFileName2 = new JLabel();
			jLabelOutFileName2.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelOutFileName2.setBounds(0, 34, 163, 25);
			panelOutFiles.add(jLabelOutFileName2);
			jLabelOutFileName2.setText("Output file name (full)");
			panelOutFiles.add(getJTextFieldFileNameAve());
			panelOutFiles.add(getJTextFieldFileNameFull());
			jPanelMain.add(getJTextFieldLy());
			jPanelMain.add(getJTextFieldLx());
		}
		return jPanelMain;
	}

	/**
	 * This method initializes jPanelMain	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getJPanelGradient() {
		if (jPanelGradient == null) {
			jLabelSigma = new JLabel();
			jLabelSigma.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelSigma.setEnabled(false);
			jLabelSigma.setBounds(new Rectangle(0, 25, 188, 20));
			jLabelSigma.setText("Sigma (gaussian std), mm");
			jLabelGradRate = new JLabel();
			jLabelGradRate.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelGradRate.setEnabled(false);
			jLabelGradRate.setBounds(new Rectangle(0, 100, 176, 20));
			jLabelGradRate.setText("Rate of attr. adding, mM/s");
			jLabelGradCenterXY = new JLabel("Source position (Xs,Ys), mm");
			jLabelGradCenterXY.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelGradCenterXY.setEnabled(false);
			jLabelGradCenterXY.setBounds(new Rectangle(0, 0, 188, 20));
			jLabelGradMin = new JLabel("Min [attractant] mM");
			jLabelGradMin.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelGradMin.setEnabled(false);
			jLabelGradMin.setDisplayedMnemonic(KeyEvent.VK_UNDEFINED);
			jLabelGradMin.setText("Min [attractant], mM");
			jLabelGradMin.setBounds(new Rectangle(0, 75, 163, 20));
			jLabelGradMin.setVisible(true);
			jLabelGauss1 = new JLabel("Max [attractant], mM");
			jLabelGauss1.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelGauss1.setEnabled(false);
			jLabelGauss1.setBounds(new Rectangle(0, 54, 163, 20));
			jLabelGauss1.setVisible(true);
			jPanelGradient=new JPanel();
			jPanelGradient.setFont(new Font("Tahoma", Font.PLAIN, 14));
			jPanelGradient.setLayout(null);
			jLabel15 = new JLabel("Gradient shape");
			jLabel15.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel15.setBounds(new Rectangle(194, 34, 105, 25));
			jLabel123 = new JLabel("t (removal), s");
			jLabel123.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel123.setEnabled(false);
			jLabel123.setBounds(new Rectangle(0, 50, 127, 20));
			jLabel123.setVisible(true);
			jLabel122 = new JLabel("t (addition), s");
			jLabel122.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel122.setEnabled(false);
			jLabel122.setBounds(new Rectangle(0, 25, 127, 20));
			jLabel122.setVisible(true);
			jLabel121 = new JLabel("Attractant step, mM");
			jLabel121.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel121.setEnabled(false);
			jLabel121.setBounds(new Rectangle(0, 0, 174, 18));
			jLabel121.setVisible(true);
			jPanelGradient.add(getJComboBoxGradient(), null);
			jPanelGradient.add(jLabel15, null);
			jPanelGradient.add(getJPanelGradParam2(), null);
			jPanelGradient.add(getJPanelGradParam1(), null);
	}
		return jPanelGradient;
	}
	/**
	 * This method initializes jPanelAdvanced	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getJPanelAdvanced() {
		if (jPanelAdvanced == null) {
			jLabelBoundary = new JLabel();
			jLabelBoundary.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelBoundary.setBounds(new Rectangle(323, 57, 146, 20));
			jLabelBoundary.setText("Boundary condition");
			jLabel11111121 = new JLabel();
			jLabel11111121.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel11111121.setBounds(new Rectangle(36, 161, 121, 20));
			jLabel11111121.setText("Resting mb(CCW)");
			jLabel1111111 = new JLabel();
			jLabel1111111.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel1111111.setBounds(new Rectangle(376, 162, 145, 20));
			jLabel1111111.setText("Cell speed (mcm/s)");
			jLabel111112 = new JLabel();
			jLabel111112.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel111112.setBounds(new Rectangle(376, 125, 146, 20));
			jLabel111112.setText("Rot. diffusion (rad2/s)");
			jLabel111111 = new JLabel();
			jLabel111111.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel111111.setBounds(new Rectangle(36, 125, 121, 20));
			jLabel111111.setText("N flagellar motors");
			jLabel16 = new JLabel();
			jLabel16.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel16.setBounds(new Rectangle(36, 54, 84, 20));
			jLabel16.setText("Medium");
			jPanelAdvanced = new JPanel();
			jPanelAdvanced.add(getJComboBoxMedium(), null);
			jPanelAdvanced.add(jLabel16, null);
			jPanelAdvanced.setLayout(null);
			jPanelAdvanced.add(jLabel111111, null);
			jPanelAdvanced.add(getJTextFieldNmotors(), null);
			jPanelAdvanced.add(jLabel111112, null);
			jPanelAdvanced.add(getJTextFieldDr(), null);
			jPanelAdvanced.add(jLabel1111111, null);
			jPanelAdvanced.add(getJTextFieldVswim(), null);
			jPanelAdvanced.add(jLabel11111121, null);
			jPanelAdvanced.add(getJTextFieldMbias0(), null);

			jPanelAdvanced.add(jLabelBoundary, null);
			jPanelAdvanced.add(getJComboBoxBoundary(), null);
		}
		return jPanelAdvanced;
	}

	/**
	 * This method initializes jTextFieldNcells	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldNcells() {
		if (jTextFieldNcells == null) {
			jTextFieldNcells = new JTextField("5");
			jTextFieldNcells.setFont(new Font("Tahoma", Font.PLAIN, 14));
			jTextFieldNcells.setBounds(new Rectangle(410, 46, 85, 20));
			jTextFieldNcells.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldNcells.setToolTipText("Start with 1 to get familiar with output");
		}
		return jTextFieldNcells;
	}

	/**
	 * This method initializes jTextFieldTmax	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldTmax() {
		if (jTextFieldTmax == null) {
			jTextFieldTmax = new JTextField("10000");
			jTextFieldTmax.setFont(new Font("Tahoma", Font.PLAIN, 14));
			jTextFieldTmax.setBounds(new Rectangle(410, 82, 85, 20));
			jTextFieldTmax.setHorizontalAlignment(JTextField.RIGHT);
		}
		return jTextFieldTmax;
	}
	/**
	 * This method initializes jTextFieldLx	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldLx() {
		if (jTextFieldLx == null) {
			jTextFieldLx = new JTextField("20");
			jTextFieldLx.setBounds(241, 226, 37, 20);
			jTextFieldLx.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldLx.setToolTipText("Arena's width, mm");
			jTextFieldLx.setHorizontalAlignment(JTextField.RIGHT);
		}
		return jTextFieldLx;
	}

	/**
	 * This method initializes jComboBoxIniPosition	
	 * 	
	 * @return javax.swing.JComboBox	
	 */
	private JComboBox getJComboBoxIniPosition() {
		if (jComboBoxIniPosition == null) {
			String[] IniPosition = {"Center of arena", "Center of left wall", "Custom (x,y)", "Random (x,y)"};
			jComboBoxIniPosition = new JComboBox(IniPosition);
			jComboBoxIniPosition.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jComboBoxIniPosition.setBounds(154, 31, 130, 20);
			jComboBoxIniPosition.setSelectedIndex(1);
			jComboBoxIniPosition.setBackground(Color.white);
			jComboBoxIniPosition.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					if(jComboBoxIniPosition.getSelectedIndex()==0){//center
						jTextFieldIniPosX.setText(String.valueOf(Double.parseDouble(jTextFieldLx.getText())/2));
						jTextFieldIniPosX.setEnabled(false);
						jTextFieldIniPosY.setText(String.valueOf(Double.parseDouble(jTextFieldLy.getText())/2));
						jTextFieldIniPosY.setEnabled(false);
					}
					else if(jComboBoxIniPosition.getSelectedIndex()==1){//left wall
						jTextFieldIniPosX.setText("0.02");
						jTextFieldIniPosX.setEnabled(false);
						jTextFieldIniPosY.setText(String.valueOf(Double.parseDouble(jTextFieldLy.getText())/2));
						jTextFieldIniPosY.setEnabled(false);
					}
					else if(jComboBoxIniPosition.getSelectedIndex()==2){//custom
						jTextFieldIniPosX.setText("");
						jTextFieldIniPosX.setEnabled(true);
						jTextFieldIniPosY.setText("");
						jTextFieldIniPosY.setEnabled(true);
					}
					else if(jComboBoxIniPosition.getSelectedIndex()==3){//random
						jTextFieldIniPosX.setText("-1");
						jTextFieldIniPosX.setEnabled(false);
						jTextFieldIniPosY.setText("-1");
						jTextFieldIniPosY.setEnabled(false);
					}
				}
			});
		}
		return jComboBoxIniPosition;
	}

	private void setJComboBoxExponentialDisplacementGradient() {
		textAreaConsole.append("\n E(x)=Rate*exp(x)");
		jPanelGradParam1.setVisible(true);
		jPanelGradParam2.setVisible(false);
		jLabelSigma.setEnabled(false);
		jTextFieldSigma.setEnabled(false);
		jLabelGradMin.setEnabled(false);
		jTextFieldMinAsp.setEnabled(false);
		jLabelGauss1.setEnabled(false);
		jTextFieldGradCenterX.setEnabled(false);
		jTextFieldGradCenterY.setEnabled(false);
		jLabelGradRate.setEnabled(true);
		jTextFieldGradRate.setEnabled(true);
	}

	/**
	 * This method initializes jComboBoxGradient	
	 * 	
	 * @return javax.swing.JComboBox	
	 */
	private JComboBox getJComboBoxGradient() {
		if (jComboBoxGradient == null) {
			String[] Gradient = {"No gradient", //0
					"Constant-activity C(x,y)", //1
					"Constant-activity C(x)", //2
					"Constant-activity C(t)", //3
					"Linear F(x)",			//4
					"Gaussian G(x,y)",		//5
					"Gaussian G(x)", 		//6
					"Exponential E(t)", 	//7
					"Exponential E(x)", 	//8
					"Stepwise add and remove"};//9
			jComboBoxGradient = new JComboBox(Gradient);
			jComboBoxGradient.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jComboBoxGradient.setBounds(new Rectangle(306, 34, 210, 25));
			jComboBoxGradient.setSelectedIndex(8);
			jComboBoxGradient.setBackground(Color.white);
			jComboBoxGradient.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					if(jComboBoxGradient.getSelectedIndex()==0){
						textAreaConsole.append("\n no attractant");
						jPanelGradParam1.setVisible(false);
						jPanelGradParam2.setVisible(false);
					}
					else if(jComboBoxGradient.getSelectedIndex()==1){//Constant-activity C(x,y)
						textAreaConsole.append("\n C(x,y)=K*C*r/((K_on-Koff)/K-C*r), r^2=(x-Lx/2)^2+(y-Ly/2)^2, K and C dep. on Lx and K_on{off}");
						jPanelGradParam1.setVisible(false);
						jPanelGradParam2.setVisible(false);
					}
					else if(jComboBoxGradient.getSelectedIndex()==2){//Constant-activity C(x)
						textAreaConsole.append("\n C(x)=K*C*x/((K_on-Koff)/K-C*x), where K and C constant, dep. on Lx and K_on{off}");
						jPanelGradParam1.setVisible(false);
						jPanelGradParam2.setVisible(false);
					}
					else if(jComboBoxGradient.getSelectedIndex()==3){//Constant-activity C(t)
						textAreaConsole.append("\n C(t)=K*C*t/((K_on-Koff)/K-C*t), where K and C constant, dep. on Lx and K_on{off}");
						jPanelGradParam1.setVisible(false);
						jPanelGradParam2.setVisible(false);
					}
					else if(jComboBoxGradient.getSelectedIndex()==4){// linear
						textAreaConsole.append("\n F(x)=(Max-Min)*x/Lx");
						jPanelGradParam1.setVisible(true);
						jPanelGradParam2.setVisible(false);
						jLabelGradMin.setEnabled(true);
						jTextFieldMinAsp.setEnabled(true);
						jLabelGauss1.setEnabled(true);
						jTextFieldMaxAsp.setEnabled(true);	
					}
					else if(jComboBoxGradient.getSelectedIndex()==5){//Gauss(x,y)
						textAreaConsole.append("\n G(x,y)=Max*exp(-r^2/(2*sigma^2), r^2=(x-Xs)^2+(y-Ys)^2");
						jPanelGradParam1.setVisible(true);
						jPanelGradParam2.setVisible(false);
						jLabelSigma.setEnabled(true);
						jTextFieldSigma.setEnabled(true);
						jLabelGradMin.setEnabled(false);
						jTextFieldMinAsp.setEnabled(false);
						jLabelGauss1.setEnabled(true);
						jTextFieldMaxAsp.setEnabled(true);
						jTextFieldGradCenterX.setEnabled(true);
						jTextFieldGradCenterY.setEnabled(true);
						jLabelGradCenterXY.setEnabled(true);
					}
					else if(jComboBoxGradient.getSelectedIndex()==6){//Gauss(x)
						textAreaConsole.append("\n G(x)=Max*exp(-(x-Xs)^2/(2*sigma^2)");
						jPanelGradParam1.setVisible(true);
						jPanelGradParam2.setVisible(false);
						jLabelSigma.setEnabled(true);
						jTextFieldSigma.setEnabled(true);
						jLabelGradMin.setEnabled(false);
						jTextFieldMinAsp.setEnabled(false);
						jLabelGauss1.setEnabled(true);
						jTextFieldMaxAsp.setEnabled(true);
						jTextFieldGradCenterX.setEnabled(true);
						jTextFieldGradCenterY.setEnabled(false);
						jLabelGradCenterXY.setEnabled(true);
					}
					else if(jComboBoxGradient.getSelectedIndex()==7){//exp(t)
						textAreaConsole.append("\n E(t)=Rate*exp(t)");
						jPanelGradParam1.setVisible(true);
						jPanelGradParam2.setVisible(false);
						jLabelSigma.setEnabled(false);
						jTextFieldSigma.setEnabled(false);
						jLabelGradMin.setEnabled(false);
						jTextFieldMinAsp.setEnabled(false);	
						jTextFieldMaxAsp.setEnabled(false);
						jLabelGradRate.setEnabled(true);	
						jTextFieldGradRate.setEnabled(true);	
						jTextFieldGradCenterX.setEnabled(false);
						jTextFieldGradCenterY.setEnabled(false);
						jLabelGradCenterXY.setEnabled(false);
						jLabelGauss1.setEnabled(false);
					}
					else if(jComboBoxGradient.getSelectedIndex()==8){//exp(x)
						setJComboBoxExponentialDisplacementGradient();
					}
					else if(jComboBoxGradient.getSelectedIndex()==9){
						textAreaConsole.append("\n Attractant stepwise added and removed");
						jPanelGradParam1.setVisible(false);
						jPanelGradParam2.setVisible(true);
						jTextFieldStep.setEnabled(true);
						jLabel121.setEnabled(true);
						jTextFieldAdditionTime.setEnabled(true);
						jLabel122.setEnabled(true);
						jTextFieldRemovalTime.setEnabled(true);
						jLabel123.setEnabled(true);
						jTextFieldGradCenterX.setEnabled(true);
						jTextFieldGradCenterY.setEnabled(true);
					}
					else{}
				}
			});
		}
		return jComboBoxGradient;
	}

	private JComboBox getJComboBoxLigandMethylation() {
		if (jComboBoxLigandMethylation == null) {
			String[] LigandMethylation = {
				"Methylation Dynamics",  // 0
				"Step Function",  // 1
				"Gaussian",  // 2
				"M = 0",  // 3
				"Exponential",  // 4
				"Linear (Two Segments)",  // 5
				"Linear (One Segment)"  // 6
			};
			jComboBoxLigandMethylation = new JComboBox(LigandMethylation);
			jComboBoxLigandMethylation.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jComboBoxLigandMethylation.setBounds(new Rectangle(370, 270, 170, 25));
			jComboBoxLigandMethylation.setSelectedIndex(0);
			jComboBoxLigandMethylation.setBackground(Color.white);
			jComboBoxLigandMethylation.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					switch (jComboBoxLigandMethylation.getSelectedIndex()) {
						case 0:
							textAreaConsole.append("\n Ligand-Methylation Relation: " + LigandMethylation[0]);
							break;
						case 1:
							textAreaConsole.append("\n Ligand-Methylation Relation: " + LigandMethylation[1]);
							break;
						case 2:
							textAreaConsole.append("\n Ligand-Methylation Relation: " + LigandMethylation[2]);
							break;
						case 3:
							textAreaConsole.append("\n Ligand-Methylation Relation: " + LigandMethylation[3]);
							break;
						case 4:
							textAreaConsole.append("\n Ligand-Methylation Relation: " + LigandMethylation[4]);
							break;
						case 5:
							textAreaConsole.append("\n Ligand-Methylation Relation: " + LigandMethylation[5]);
							break;
						case 6:
							textAreaConsole.append("\n Ligand-Methylation Relation: " + LigandMethylation[6]);
							break;
					}
				}
			});
		}
		return jComboBoxLigandMethylation;
	}

	/**
	 * This method initializes jComboBoxMedium	
	 * 	
	 * @return javax.swing.JComboBox	
	 */
	private JComboBox getJComboBoxMedium() {
		if (jComboBoxMedium == null) {
			String[] Medium = {"Liquid","Agar"};
			jComboBoxMedium = new JComboBox(Medium);
			jComboBoxMedium.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jComboBoxMedium.setBounds(new Rectangle(130, 54, 110, 25));
			jComboBoxMedium.setBackground(Color.white);
			jComboBoxMedium.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					if(jComboBoxMedium.getSelectedIndex()==0){
						textAreaConsole.setBackground(Color.white);
						textAreaConsole.append("\n Simple run-tumble motion, no tricks.");
					}
					if(jComboBoxMedium.getSelectedIndex()==1){
						textAreaConsole.setBackground(Color.white);
						textAreaConsole.append("\n Motion \"run-trap-tumble\" for agar medium. Traps may occur once in ~2 s due to agar filaments.");
					}
				}
			});
		}
		return jComboBoxMedium;
	}

	/**
	 * This method initializes jTextFieldK_offTar	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldK_offTar() {
		if (jTextFieldK_offTar == null) {
			jTextFieldK_offTar = new JTextField();
			jTextFieldK_offTar.setHorizontalAlignment(SwingConstants.RIGHT);
			jTextFieldK_offTar.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldK_offTar.setToolTipText("Diss. constant K_off(Tar)");
			jTextFieldK_offTar.setBounds(new Rectangle(58, 40, 61, 20));
			jTextFieldK_offTar.setText("0.02");
		}
		return jTextFieldK_offTar;
	}

	/**
	 * This method initializes jTextFieldK_onTar	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldK_onTar() {
		if (jTextFieldK_onTar == null) {
			jTextFieldK_onTar = new JTextField();
			jTextFieldK_onTar.setHorizontalAlignment(SwingConstants.RIGHT);
			jTextFieldK_onTar.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldK_onTar.setToolTipText("Diss. constant K_on(Tar)");
			jTextFieldK_onTar.setBounds(new Rectangle(58, 68, 61, 20));
			jTextFieldK_onTar.setText("0.5");
		}
		return jTextFieldK_onTar;
	}

	/**
	 * This method initializes jTextFieldK_offTsr	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldK_offTsr() {
		if (jTextFieldK_offTsr == null) {
			jTextFieldK_offTsr = new JTextField();
			jTextFieldK_offTsr.setHorizontalAlignment(SwingConstants.RIGHT);
			jTextFieldK_offTsr.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldK_offTsr.setToolTipText("Diss. constant K_off(Tsr)");
			jTextFieldK_offTsr.setBounds(new Rectangle(163, 40, 76, 20));
			jTextFieldK_offTsr.setText("100");
		}
		return jTextFieldK_offTsr;
	}

	/**
	 * This method initializes jTextFieldK_onTsr	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldK_onTsr() {
		if (jTextFieldK_onTsr == null) {
			jTextFieldK_onTsr = new JTextField();
			jTextFieldK_onTsr.setHorizontalAlignment(SwingConstants.RIGHT);
			jTextFieldK_onTsr.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldK_onTsr.setToolTipText("Diss. const. K_on(Tsr)");
			jTextFieldK_onTsr.setBounds(new Rectangle(163, 68, 76, 20));
			jTextFieldK_onTsr.setText("1000000");
		}
		return jTextFieldK_onTsr;
	}

	/**
	 * This method initializes jComboBoxModelRun	
	 * 	
	 * @return javax.swing.JComboBox	
	 */
	private JComboBox getJComboBoxModelRun() {
		if (jComboBoxModelRun == null) {
			String[] ModelRun = {"Voting motors","Distortion factor (?)"};
			jComboBoxModelRun = new JComboBox(ModelRun);
			jComboBoxModelRun.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jComboBoxModelRun.setBounds(new Rectangle(195, 61, 250, 25));
			jComboBoxModelRun.setSelectedIndex(0);
			jComboBoxModelRun.setBackground(Color.white);
			jComboBoxModelRun.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					if(jComboBoxModelRun.getSelectedIndex()==0){
						textAreaConsole.setBackground(Color.white);
						textAreaConsole.append("\n Run->Tumble when majority of motors become CW. Tumble->Run when majority CCW.");
					}
					if(jComboBoxModelRun.getSelectedIndex()==1){
						textAreaConsole.setBackground(Color.white);
						textAreaConsole.append("\n More complicated, see (Vladimirov, Lebiedz, Sourjik, 2010).");
					}
				}
			});
		}
		return jComboBoxModelRun;
	}

	/**
	 * This method initializes jComboBoxModelTumble	
	 * 	
	 * @return javax.swing.JComboBox	
	 */
	private JComboBox getJComboBoxModelTumble() {
		if (jComboBoxModelTumble == null) {
			String[] ModelTumble = {"Isotropic","Anisotropic (?)"};
			jComboBoxModelTumble = new JComboBox(ModelTumble);
			jComboBoxModelTumble.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jComboBoxModelTumble.setBounds(new Rectangle(195, 111, 250, 25));
			jComboBoxModelTumble.setBackground(Color.white);
			jComboBoxModelTumble.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					if(jComboBoxModelTumble.getSelectedIndex()==0){
						textAreaConsole.setBackground(Color.white);
						textAreaConsole.append("\n Tumbling angle is random with p.d.f. p(a)=0.5(1+Cos(a))Sin(a), average <a>=67.5 deg.");
					}
					if(jComboBoxModelTumble.getSelectedIndex()==1){
						textAreaConsole.setBackground(Color.white);
						textAreaConsole.append("\n Angle is random and depends on N(CW motors), <a>=67.5 deg. (Vladimirov et al., 2010)");
					}
				}
			});
		}
		return jComboBoxModelTumble;
	}

	/**
	 * This method initializes jTextFieldAdRate	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldAdRate() {
		if (jTextFieldAdRate == null) {
			jTextFieldAdRate = new JTextField();
			jTextFieldAdRate.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldAdRate.setBounds(new Rectangle(438, 96, 59, 20));
			jTextFieldAdRate.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldAdRate.setText("1.0");
		}
		return jTextFieldAdRate;
	}

	/**
	 * This method initializes jTextFieldNmotors	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldNmotors() {
		if (jTextFieldNmotors == null) {
			jTextFieldNmotors = new JTextField();
			jTextFieldNmotors.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldNmotors.setBounds(new Rectangle(180, 126, 60, 20));
			jTextFieldNmotors.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldNmotors.setToolTipText("Number of flagellar motors per cell");
			jTextFieldNmotors.setText("5");
		}
		return jTextFieldNmotors;
	}

	/**
	 * This method initializes jTextFieldDr	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldDr() {
		if (jTextFieldDr == null) {
			jTextFieldDr = new JTextField();
			jTextFieldDr.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldDr.setBounds(new Rectangle(534, 124, 56, 20));
			jTextFieldDr.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldDr.setToolTipText("Rotational diffusion due to Brownian motion. Depends on medium temperature and viscosity.");
			jTextFieldDr.setText("0.062");
		}
		return jTextFieldDr;
	}

	/**
	 * This method initializes jTextFieldVswim	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldVswim() {
		if (jTextFieldVswim == null) {
			jTextFieldVswim = new JTextField();
			jTextFieldVswim.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldVswim.setBounds(new Rectangle(553, 161, 37, 20));
			jTextFieldVswim.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldVswim.setText("20");
		}
		return jTextFieldVswim;
	}

	/**
	 * This method initializes jTextFieldNtar	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldNtar() {
		if (jTextFieldNtar == null) {
			jTextFieldNtar = new JTextField();
			jTextFieldNtar.setBounds(0, 48, 37, 20);
			jTextFieldNtar.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldNtar.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldNtar.setToolTipText("Tar receptors in cluster");
			jTextFieldNtar.setText("6");
		}
		return jTextFieldNtar;
	}

	/**
	 * This method initializes jTextFieldNtsr	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldNtsr() {
		if (jTextFieldNtsr == null) {
			jTextFieldNtsr = new JTextField();
			jTextFieldNtsr.setBounds(85, 48, 37, 20);
			jTextFieldNtsr.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldNtsr.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldNtsr.setToolTipText("Tsr receptors in cluster");
			jTextFieldNtsr.setText("12");
		}
		return jTextFieldNtsr;
	}

	/**
	 * This method initializes jTextFieldMbias0	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldMbias0() {
		if (jTextFieldMbias0 == null) {
			jTextFieldMbias0 = new JTextField();
			jTextFieldMbias0.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldMbias0.setBounds(new Rectangle(180, 161, 60, 20));
			jTextFieldMbias0.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldMbias0.setToolTipText("CCW motor bias in adapted state. Values between 0 and 1.");
			jTextFieldMbias0.setText("0.65");
		}
		return jTextFieldMbias0;
	}

	/**
	 * This method initializes jProgressBar	
	 * 	
	 * @return javax.swing.JProgressBar	
	 */
	private JProgressBar getJProgressBar() {
		if (jProgressBar == null) {
			jProgressBar = new JProgressBar();
			jProgressBar.setPreferredSize(new Dimension(760, 16));
			jProgressBar.setBounds(new Rectangle(12, 532, 760, 16));
		}
		return jProgressBar;
	}

	/**
	 * This method initializes jPanel	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getJPanelNetwork() {
		if (jPanelNetwork == null) {
			jLabel117 = new JLabel();
			jLabel117.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel117.setBounds(new Rectangle(65, 239, 166, 20));
			jLabel117.setText("[CheZ] relative to wt, %");
			jLabel116 = new JLabel();
			jLabel116.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel116.setBounds(new Rectangle(65, 209, 166, 20));
			jLabel116.setText("[CheY] relative to wt, %");
			jLabel1141 = new JLabel();
			jLabel1141.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel1141.setBounds(new Rectangle(65, 159, 166, 20));
			jLabel1141.setText("[CheB] relative to wt, %");
			jLabel114 = new JLabel();
			jLabel114.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel114.setBounds(new Rectangle(65, 129, 166, 20));
			jLabel114.setText("[CheR] relative to wt, %");
			jLabel113 = new JLabel();
			jLabel113.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel113.setText("Adaptation precision, %");
			jLabel113.setBounds(new Rectangle(65, 96, 166, 20));
			jPanelNetwork = new JPanel();
			jPanelNetwork.setLayout(null);
			jLabel11111 = new JLabel();
			jLabel11111.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabel11111.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel11111.setBounds(new Rectangle(289, 96, 137, 20));
			jLabel11111.setText("Adapt. rate, relative");
			jPanelNetwork.add(jLabel11111, null);
			jPanelNetwork.add(getJTextFieldAdRate(), null);
			jPanelNetwork.add(jLabel113, null);
			jPanelNetwork.add(getJTextFieldAdPrecision(), null);
			jPanelNetwork.add(getJpanelCluster());
			jPanelNetwork.add(jLabel114, null);
			jPanelNetwork.add(jLabel1141, null);
			jPanelNetwork.add(getJTextFieldCheR(), null);
			jPanelNetwork.add(getJTextFieldCheB(), null);
			jPanelNetwork.add(getJTextFieldCheY(), null);
			jPanelNetwork.add(getJTextFieldCheZ(), null);
			jPanelNetwork.add(jLabel116, null);
			jPanelNetwork.add(jLabel117, null);
			jPanelNetwork.add(getJPanelKd(), null);
			
			JLabel lblcheaRelativeTo = new JLabel("[CheA] relative to wt, %");
			lblcheaRelativeTo.setFont(new Font("Tahoma", Font.PLAIN, 15));
			lblcheaRelativeTo.setBounds(65, 47, 166, 20);
			jPanelNetwork.add(lblcheaRelativeTo);
			
			jTextFieldCheA = new JTextField();
			jTextFieldCheA.setHorizontalAlignment(SwingConstants.RIGHT);
			jTextFieldCheA.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldCheA.setText("100");
			jTextFieldCheA.setBounds(240, 47, 46, 20);
			jPanelNetwork.add(jTextFieldCheA);
			jTextFieldCheA.setColumns(10);
			
			textFieldkA = new JTextField();
			textFieldkA.setHorizontalAlignment(SwingConstants.RIGHT);
			textFieldkA.setFont(new Font("Tahoma", Font.PLAIN, 15));
			textFieldkA.setText("5");
			textFieldkA.setBounds(438, 47, 59, 20);
			jPanelNetwork.add(textFieldkA);
			textFieldkA.setColumns(10);
			
			JLabel lblAutophosphRate = new JLabel("Autophosph. rate");
			lblAutophosphRate.setHorizontalAlignment(SwingConstants.RIGHT);
			lblAutophosphRate.setFont(new Font("Tahoma", Font.PLAIN, 15));
			lblAutophosphRate.setBounds(306, 46, 120, 22);
			jPanelNetwork.add(lblAutophosphRate);
			
			textFieldkY = new JTextField();
			textFieldkY.setHorizontalAlignment(SwingConstants.RIGHT);
			textFieldkY.setFont(new Font("Tahoma", Font.PLAIN, 15));
			textFieldkY.setText("100");
			textFieldkY.setBounds(438, 209, 59, 20);
			jPanelNetwork.add(textFieldkY);
			textFieldkY.setColumns(10);
			
			textFieldkZ = new JTextField();
			textFieldkZ.setFont(new Font("Tahoma", Font.PLAIN, 15));
			textFieldkZ.setHorizontalAlignment(SwingConstants.RIGHT);
			textFieldkZ.setText("30");
			textFieldkZ.setBounds(438, 239, 59, 20);
			jPanelNetwork.add(textFieldkZ);
			textFieldkZ.setColumns(10);
			
			JLabel lblCheyPhosphRate = new JLabel("CheY phosph. rate");
			lblCheyPhosphRate.setHorizontalAlignment(SwingConstants.RIGHT);
			lblCheyPhosphRate.setFont(new Font("Tahoma", Font.PLAIN, 15));
			lblCheyPhosphRate.setBounds(299, 209, 127, 20);
			jPanelNetwork.add(lblCheyPhosphRate);
			
			JLabel lblChezPhosphRate = new JLabel("CheZ phosph. rate");
			lblChezPhosphRate.setHorizontalAlignment(SwingConstants.RIGHT);
			lblChezPhosphRate.setFont(new Font("Tahoma", Font.PLAIN, 15));
			lblChezPhosphRate.setBounds(298, 239, 127, 20);
			jPanelNetwork.add(lblChezPhosphRate);
			
			textFieldkR = new JTextField();
			textFieldkR.setFont(new Font("Tahoma", Font.PLAIN, 15));
			textFieldkR.setHorizontalAlignment(SwingConstants.RIGHT);
			textFieldkR.setText("0.0182");
			textFieldkR.setBounds(438, 129, 59, 20);
			jPanelNetwork.add(textFieldkR);
			textFieldkR.setColumns(10);
			
			JLabel lblCherCatalyticRate = new JLabel("CheR catalytic rate");
			lblCherCatalyticRate.setHorizontalAlignment(SwingConstants.RIGHT);
			lblCherCatalyticRate.setFont(new Font("Tahoma", Font.PLAIN, 15));
			lblCherCatalyticRate.setBounds(298, 129, 128, 20);
			jPanelNetwork.add(lblCherCatalyticRate);
			
			textFieldkB = new JTextField();
			textFieldkB.setFont(new Font("Tahoma", Font.PLAIN, 15));
			textFieldkB.setHorizontalAlignment(SwingConstants.RIGHT);
			textFieldkB.setText("0.0364");
			textFieldkB.setBounds(438, 159, 59, 20);
			jPanelNetwork.add(textFieldkB);
			textFieldkB.setColumns(10);
			
			JLabel lblChebCatalyticRate = new JLabel("CheB catalytic rate");
			lblChebCatalyticRate.setHorizontalAlignment(SwingConstants.RIGHT);
			lblChebCatalyticRate.setFont(new Font("Tahoma", Font.PLAIN, 15));
			lblChebCatalyticRate.setBounds(298, 159, 128, 20);
			jPanelNetwork.add(lblChebCatalyticRate);

			JLabel jLabelLigandMethylation = new JLabel("Ligand-Methylation Relation");
			jLabelLigandMethylation.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelLigandMethylation.setBounds(new Rectangle(175, 270, 200, 25));
			jPanelNetwork.add(jLabelLigandMethylation);

			jPanelNetwork.add(getJComboBoxLigandMethylation(), null);
		}
		return jPanelNetwork;
	}

	/**
	 * This method initializes jTextFieldAdPrecision	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldAdPrecision() {
		if (jTextFieldAdPrecision == null) {
			jTextFieldAdPrecision = new JTextField();
			jTextFieldAdPrecision.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldAdPrecision.setBounds(new Rectangle(240, 96, 46, 20));
			jTextFieldAdPrecision.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldAdPrecision.setText("100");
		}
		return jTextFieldAdPrecision;
	}

	/**
	 * This method initializes jTextFieldCheR	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldCheR() {
		if (jTextFieldCheR == null) {
			jTextFieldCheR = new JTextField();
			jTextFieldCheR.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldCheR.setBounds(new Rectangle(240, 129, 46, 20));
			jTextFieldCheR.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldCheR.setText("100");
		}
		return jTextFieldCheR;
	}

	/**
	 * This method initializes jTextFieldCheB	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldCheB() {
		if (jTextFieldCheB == null) {
			jTextFieldCheB = new JTextField();
			jTextFieldCheB.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldCheB.setBounds(new Rectangle(240, 159, 46, 20));
			jTextFieldCheB.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldCheB.setText("100");
		}
		return jTextFieldCheB;
	}

	/**
	 * This method initializes jComboBoxModelMotor	
	 * 	
	 * @return javax.swing.JComboBox	
	 */
	private JComboBox getJComboBoxModelMotor() {
		if (jComboBoxModelMotor == null) {
			String[] ModelMotor = {"T_ccw(CheYp), T_cw is constant","T_ccw(CheYp) and T_cw(CheYp) (?)"};
			jComboBoxModelMotor = new JComboBox(ModelMotor);
			jComboBoxModelMotor.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jComboBoxModelMotor.setBounds(new Rectangle(195, 161, 250, 25));
			jComboBoxModelMotor.setSelectedIndex(0);
			jComboBoxModelMotor.setBackground(Color.white);
			jComboBoxModelMotor.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					if(jComboBoxModelMotor.getSelectedIndex()==0){
						textAreaConsole.setBackground(Color.white);
						textAreaConsole.append("\n Mean CCW interval depends on CheY-P, but mean CW interval does not.");
					}
					if(jComboBoxModelMotor.getSelectedIndex()==1){
						textAreaConsole.setBackground(Color.white);
						textAreaConsole.append("\n Mean CCW and CW intervals depend on CheY-P. See Vladimirov, PhD Thesis, Appendix B");
					}
				}
			});
		}
		return jComboBoxModelMotor;
	}

	/**
	 * This method initializes jTextFieldCheY	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldCheY() {
		if (jTextFieldCheY == null) {
			jTextFieldCheY = new JTextField();
			jTextFieldCheY.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldCheY.setBounds(new Rectangle(240, 209, 46, 20));
			jTextFieldCheY.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldCheY.setText("100");
		}
		return jTextFieldCheY;
	}

	/**
	 * This method initializes jTextFieldCheZ	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldCheZ() {
		if (jTextFieldCheZ == null) {
			jTextFieldCheZ = new JTextField();
			jTextFieldCheZ.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldCheZ.setBounds(new Rectangle(240, 239, 46, 20));
			jTextFieldCheZ.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldCheZ.setText("100");
		}
		return jTextFieldCheZ;
	}

	/**
	 * This method initializes jTextFieldMeth	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldMeth() {
		if (jTextFieldMeth == null) {
			jTextFieldMeth = new JTextField();
			jTextFieldMeth.setBounds(120, 81, 42, 20);
			jTextFieldMeth.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldMeth.setHorizontalAlignment(SwingConstants.RIGHT);
			jTextFieldMeth.setText("1.924");
		}
		return jTextFieldMeth;
	}

	/**
	 * This method initializes jTextFieldStep	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldStep() {
		if (jTextFieldStep == null) {
			jTextFieldStep = new JTextField();
			jTextFieldStep.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldStep.setText("0.03");
			jTextFieldStep.setToolTipText("Step of attractant concentration, Asp for example");
			jTextFieldStep.setEditable(true);
			jTextFieldStep.setEnabled(false);
			jTextFieldStep.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldStep.setBounds(new Rectangle(186, 0, 46, 20));
			jTextFieldStep.setVisible(true);
		}
		return jTextFieldStep;
	}

	/**
	 * This method initializes jTextFieldAdditionTime	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldAdditionTime() {
		if (jTextFieldAdditionTime == null) {
			jTextFieldAdditionTime = new JTextField();
			jTextFieldAdditionTime.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldAdditionTime.setText("50");
			jTextFieldAdditionTime.setEnabled(false);
			jTextFieldAdditionTime.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldAdditionTime.setBounds(new Rectangle(186, 25, 46, 20));
			jTextFieldAdditionTime.setVisible(true);
		}
		return jTextFieldAdditionTime;
	}

	/**
	 * This method initializes jTextFieldRemovalTime	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldRemovalTime() {
		if (jTextFieldRemovalTime == null) {
			jTextFieldRemovalTime = new JTextField();
			jTextFieldRemovalTime.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldRemovalTime.setText("350");
			jTextFieldRemovalTime.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldRemovalTime.setEnabled(false);
			jTextFieldRemovalTime.setBounds(new Rectangle(186, 50, 46, 20));
			jTextFieldRemovalTime.setVisible(true);
		}
		return jTextFieldRemovalTime;
	}

	/**
	 * This method initializes jTextFieldMaxAsp	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldMaxAsp() {
		if (jTextFieldMaxAsp == null) {
			jTextFieldMaxAsp = new JTextField("0.1");
			jTextFieldMaxAsp.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldMaxAsp.setToolTipText("Set MAX concentration of attractant, [mM]");
			jTextFieldMaxAsp.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldMaxAsp.setEnabled(false);
			jTextFieldMaxAsp.setBounds(new Rectangle(199, 50, 37, 20));
			jTextFieldMaxAsp.setVisible(true);
		}
		return jTextFieldMaxAsp;
	}

	/**
	 * This method initializes jTextFieldMinAsp	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldMinAsp() {
		if (jTextFieldMinAsp == null) {
			jTextFieldMinAsp = new JTextField("0");
			jTextFieldMinAsp.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldMinAsp.setToolTipText("Set MIN concentration of attractant, [mM]");
			jTextFieldMinAsp.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldMinAsp.setEditable(true);
			jTextFieldMinAsp.setEnabled(false);
			jTextFieldMinAsp.setBounds(new Rectangle(199, 75, 37, 20));
			jTextFieldMinAsp.setVisible(true);
		}
		return jTextFieldMinAsp;
	}

	/**
	 * This method initializes jTextFieldLy	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldLy() {
		if (jTextFieldLy == null) {
			jTextFieldLy = new JTextField();
			jTextFieldLy.setBounds(35, 19, 37, 20);
			jTextFieldLy.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldLy.setToolTipText("Arena's height, mm");
			jTextFieldLy.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldLy.setText("20");
		}
		return jTextFieldLy;
	}

	/**
	 * This method initializes jTextFieldGradCenterX	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldGradCenterX() {
		if (jTextFieldGradCenterX == null) {
			jTextFieldGradCenterX = new JTextField();
			jTextFieldGradCenterX.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldGradCenterX.setText("10");
			jTextFieldGradCenterX.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldGradCenterX.setEditable(true);
			jTextFieldGradCenterX.setEnabled(false);
			jTextFieldGradCenterX.setBounds(new Rectangle(199, 0, 37, 20));
			jTextFieldGradCenterX.setToolTipText("X position of the center of gradient");
		}
		return jTextFieldGradCenterX;
	}

	/**
	 * This method initializes jTextFieldGradCenterY	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldGradCenterY() {
		if (jTextFieldGradCenterY == null) {
			jTextFieldGradCenterY = new JTextField();
			jTextFieldGradCenterY.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldGradCenterY.setText("10");
			jTextFieldGradCenterY.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldGradCenterY.setEnabled(false);
			jTextFieldGradCenterY.setBounds(new Rectangle(245, 0, 37, 20));
			jTextFieldGradCenterY.setToolTipText("Y position of the gradient center");
		}
		return jTextFieldGradCenterY;
	}

	/**
	 * This method initializes jTextFieldGradRate	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldGradRate() {
		if (jTextFieldGradRate == null) {
			jTextFieldGradRate = new JTextField();
			jTextFieldGradRate.setHorizontalAlignment(SwingConstants.RIGHT);
			jTextFieldGradRate.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldGradRate.setEnabled(false);
			jTextFieldGradRate.setBounds(new Rectangle(188, 100, 48, 20));
			jTextFieldGradRate.setText("0.001");
		}
		return jTextFieldGradRate;
	}

	/**
	 * This method initializes jTextFieldSigma	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldSigma() {
		if (jTextFieldSigma == null) {
			jTextFieldSigma = new JTextField("5");
			jTextFieldSigma.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldSigma.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldSigma.setBounds(new Rectangle(199, 25, 37, 20));
			jTextFieldSigma.setEnabled(false);
		}
		return jTextFieldSigma;
	}

	/**
	 * This method initializes jPanelGradParam2	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getJPanelGradParam2() {
		if (jPanelGradParam2 == null) {
			jPanelGradParam2 = new JPanel();
			jPanelGradParam2.setLayout(null);
			jPanelGradParam2.setBounds(new Rectangle(482, 126, 244, 136));
			jPanelGradParam2.add(jLabel121, null);
			jPanelGradParam2.add(getJTextFieldStep(), null);
			jPanelGradParam2.add(getJTextFieldAdditionTime(), null);
			jPanelGradParam2.add(jLabel123, null);
			jPanelGradParam2.add(getJTextFieldRemovalTime(), null);
			jPanelGradParam2.add(jLabel122, null);
		}
		return jPanelGradParam2;
	}

	/**
	 * This method initializes jPanelGradParam1	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getJPanelGradParam1() {
		if (jPanelGradParam1 == null) {
			jPanelGradParam1 = new JPanel();
			jPanelGradParam1.setLayout(null);
			jPanelGradParam1.setBounds(new Rectangle(64, 126, 294, 136));
			jPanelGradParam1.add(jLabelGradCenterXY, null);
			jPanelGradParam1.add(getJTextFieldGradCenterX(), null);
			jPanelGradParam1.add(getJTextFieldGradCenterY(), null);
			jPanelGradParam1.add(jLabelGradRate, null);
			jPanelGradParam1.add(jLabelGauss1, null);
			jPanelGradParam1.add(jLabelGradMin, null);
			jPanelGradParam1.add(getJTextFieldGradRate(), null);
			jPanelGradParam1.add(jLabelSigma, null);
			jPanelGradParam1.add(getJTextFieldMaxAsp(), null);
			jPanelGradParam1.add(getJTextFieldMinAsp(), null);
			jPanelGradParam1.add(getJTextFieldSigma(), null);
		}
		return jPanelGradParam1;
	}

	/**
	 * This method initializes jPanelKd	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getJPanelKd() {
		if (jPanelKd == null) {
			jPanelKd = new JPanel();
			jPanelKd.setLayout(null);
			jPanelKd.setBounds(new Rectangle(509, 42, 294, 104));
			jPanelKd.add(getJTextFieldK_offTsr(), null);
			jPanelKd.add(getJTextFieldK_onTar(), null);
			jPanelKd.add(getJTextFieldK_onTsr(), null);
			jPanelKd.add(jLabel1102, null);
			jPanelKd.add(jLabel17, null);
			jPanelKd.add(getJTextFieldK_offTar(), null);
			jPanelKd.add(jLabel11011, null);
			jPanelKd.add(jLabel19, null);
			jPanelKd.add(jLabel18, null);
			jPanelKd.add(jLabel110, null);
			jPanelKd.add(jLabel1101, null);
		}
		return jPanelKd;
	}

	/**
	 * This method initializes jTextFieldIniPosX	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldIniPosX() {
		if (jTextFieldIniPosX == null) {
			jTextFieldIniPosX = new JTextField("0.02");
			jTextFieldIniPosX.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldIniPosX.setBounds(198, 0, 37, 20);
			jTextFieldIniPosX.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldIniPosX.setEnabled(false);
		}
		return jTextFieldIniPosX;
	}

	/**
	 * This method initializes jTextFieldIniPosY	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldIniPosY() {
		if (jTextFieldIniPosY == null) {
			jTextFieldIniPosY = new JTextField("10.0");
			jTextFieldIniPosY.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldIniPosY.setBounds(247, 0, 37, 20);
			jTextFieldIniPosY.setEnabled(false);
			jTextFieldIniPosY.setHorizontalAlignment(JTextField.RIGHT);
		}
		return jTextFieldIniPosY;
	}

	/**
	 * This method initializes jComboBoxIniOri	
	 * 	
	 * @return javax.swing.JComboBox	
	 */
	private JComboBox getJComboBoxIniOri() {
		if (jComboBoxIniOri == null) {
			String[] IniPosition = {"-->", "up", "<--", "down", "Random"};
			jComboBoxIniOri = new JComboBox(IniPosition);
			jComboBoxIniOri.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jComboBoxIniOri.setBounds(409, 31, 86, 20);
			jComboBoxIniOri.setSelectedIndex(4);
			jComboBoxIniOri.setBackground(Color.white);
		}
		return jComboBoxIniOri;
	}

	/**
	 * This method initializes jPanelOutput	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getJPanelOutput() {
		if (jPanelOutput == null) {
			jLabelOutDt = new JLabel();
			jLabelOutDt.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelOutDt.setText("Output time step, s");
			jLabelOutDt.setBounds(new Rectangle(548, 34, 144, 25));
			jPanelOutput = new JPanel();
			jPanelOutput.setLayout(null);
			jPanelOutput.add(jLabelOutDt, null);
			jPanelOutput.add(getJTextFieldOutDt(), null);
			jPanelOutput.add(getPanelCheckBoxes());
		}
		return jPanelOutput;
	}

	/**
	 * This method initializes jTextFieldOutDt	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldOutDt() {
		if (jTextFieldOutDt == null) {
			jTextFieldOutDt = new JTextField();
			jTextFieldOutDt.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldOutDt.setBounds(new Rectangle(704, 32, 37, 25));
			jTextFieldOutDt.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldOutDt.setText("0.1");
		}
		return jTextFieldOutDt;
	}

	/**
	 * This method initializes jTextFieldFileNameAve	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldFileNameAve() {
		if (jTextFieldFileNameAve == null) {
			jTextFieldFileNameAve = new JTextField();
			jTextFieldFileNameAve.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldFileNameAve.setBounds(177, 1, 237, 25);
			jTextFieldFileNameAve.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldFileNameAve.setText("averages.txt");
		}
		return jTextFieldFileNameAve;
	}

	/**
	 * This method initializes jTextFieldFileNameFull	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldFileNameFull() {
		if (jTextFieldFileNameFull == null) {
			jTextFieldFileNameFull = new JTextField();
			jTextFieldFileNameFull.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldFileNameFull.setBounds(177, 35, 237, 25);
			jTextFieldFileNameFull.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldFileNameFull.setText("individuals.txt");
		}
		return jTextFieldFileNameFull;
	}

	/**
	 * This method initializes jCheckBoxOutXY	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getJCheckBoxOutXY() {
		if (jCheckBoxOutXY == null) {
			jCheckBoxOutXY = new JCheckBox();
			jCheckBoxOutXY.setBounds(155, 34, 25, 25);
			jCheckBoxOutXY.setSelected(true);
		}
		return jCheckBoxOutXY;
	}

	/**
	 * This method initializes jCheckBoxOutOri	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getJCheckBoxOutOri() {
		if (jCheckBoxOutOri == null) {
			jCheckBoxOutOri = new JCheckBox();
			jCheckBoxOutOri.setBounds(155, 69, 25, 25);
			jCheckBoxOutOri.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jCheckBoxOutOri.setSelected(true);
		}
		return jCheckBoxOutOri;
	}

	/**
	 * This method initializes jCheckBoxOutCheA	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getJCheckBoxOutCheA() {
		if (jCheckBoxOutCheA == null) {
			jCheckBoxOutCheA = new JCheckBox();
			jCheckBoxOutCheA.setBounds(155, 106, 25, 25);
			jCheckBoxOutCheA.setSelected(true);
		}
		return jCheckBoxOutCheA;
	}

	/**
	 * This method initializes jCheckBoxOutCheY	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getJCheckBoxOutCheY() {
		if (jCheckBoxOutCheY == null) {
			jCheckBoxOutCheY = new JCheckBox();
			jCheckBoxOutCheY.setBounds(155, 143, 25, 25);
			jCheckBoxOutCheY.setSelected(true);
		}
		return jCheckBoxOutCheY;
	}

	/**
	 * This method initializes jCheckBoxOutMeth	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getJCheckBoxOutMeth() {
		if (jCheckBoxOutMeth == null) {
			jCheckBoxOutMeth = new JCheckBox();
			jCheckBoxOutMeth.setBounds(155, 178, 25, 25);
			jCheckBoxOutMeth.setSelected(true);
		}
		return jCheckBoxOutMeth;
	}

	/**
	 * This method initializes jCheckBoxOutMb	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getJCheckBoxOutMb() {
		if (jCheckBoxOutMb == null) {
			jCheckBoxOutMb = new JCheckBox();
			jCheckBoxOutMb.setBounds(381, 34, 25, 25);
			jCheckBoxOutMb.setSelected(true);
		}
		return jCheckBoxOutMb;
	}

	/**
	 * This method initializes jCheckBoxOutRun	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getJCheckBoxOutRun() {
		if (jCheckBoxOutRun == null) {
			jCheckBoxOutRun = new JCheckBox();
			jCheckBoxOutRun.setBounds(381, 143, 25, 25);
			jCheckBoxOutRun.setFont(new Font("Tahoma", Font.PLAIN, 15));
		}
		return jCheckBoxOutRun;
	}

	/**
	 * This method initializes jCheckBoxOutAttr	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getJCheckBoxOutAttr() {
		if (jCheckBoxOutAttr == null) {
			jCheckBoxOutAttr = new JCheckBox();
			jCheckBoxOutAttr.setBounds(381, 178, 25, 25);
			jCheckBoxOutAttr.setFont(new Font("Tahoma", Font.PLAIN, 15));
		}
		return jCheckBoxOutAttr;
	}

	/**
	 * This method initializes jCheckBoxOutNmot	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getJCheckBoxOutNmot() {
		if (jCheckBoxOutNmot == null) {
			jCheckBoxOutNmot = new JCheckBox();
			jCheckBoxOutNmot.setBounds(381, 106, 25, 25);
		}
		return jCheckBoxOutNmot;
	}

	/**
	 * This method initializes jCheckBoxOutSwitchFreq	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getJCheckBoxOutSwitchFreq() {
		if (jCheckBoxOutSwitchFreq == null) {
			jCheckBoxOutSwitchFreq = new JCheckBox();
			jCheckBoxOutSwitchFreq.setBounds(381, 69, 25, 25);
		}
		return jCheckBoxOutSwitchFreq;
	}

	/**
	 * This method initializes jCheckBoxCheckAll	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getJCheckBoxCheckAll() {
		if (jCheckBoxCheckAll == null) {
			jCheckBoxCheckAll = new JCheckBox();
			jCheckBoxCheckAll.setBounds(543, 143, 25, 25);
			jCheckBoxCheckAll.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					if(jCheckBoxCheckAll.isSelected()){
						jCheckBoxOutXY.setSelected(true);
						jCheckBoxOutOri.setSelected(true);
						jCheckBoxOutCheA.setSelected(true);
						jCheckBoxOutCheY.setSelected(true);
						jCheckBoxOutMeth.setSelected(true);
						jCheckBoxOutMb.setSelected(true);
						jCheckBoxOutSwitchFreq.setSelected(true);
						jCheckBoxOutNmot.setSelected(true);
						jCheckBoxOutRun.setSelected(true);
						jCheckBoxOutAttr.setSelected(true);
						jCheckBoxOutPon.setSelected(true);
						jCheckBoxUncheckAll.setSelected(false);
					}
				}
			});
		}
		return jCheckBoxCheckAll;
	}

	/**
	 * This method initializes jCheckBoxUncheckAll	
	 * 	
	 * @return javax.swing.JCheckBox	
	 */
	private JCheckBox getJCheckBoxUncheckAll() {
		if (jCheckBoxUncheckAll == null) {
			jCheckBoxUncheckAll = new JCheckBox();
			jCheckBoxUncheckAll.setBounds(543, 178, 25, 25);
			jCheckBoxUncheckAll.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					if(jCheckBoxUncheckAll.isSelected()){
						jCheckBoxOutXY.setSelected(false);
						jCheckBoxOutOri.setSelected(false);
						jCheckBoxOutCheA.setSelected(false);
						jCheckBoxOutCheY.setSelected(false);
						jCheckBoxOutMeth.setSelected(false);
						jCheckBoxOutMb.setSelected(false);
						jCheckBoxOutSwitchFreq.setSelected(false);
						jCheckBoxOutNmot.setSelected(false);
						jCheckBoxOutRun.setSelected(false);
						jCheckBoxOutAttr.setSelected(false);
						jCheckBoxOutPon.setSelected(false);
						jCheckBoxCheckAll.setSelected(false);
					}
				}
			});
		}
		return jCheckBoxUncheckAll;
	}

	/**
	 * This method initializes jPanelAlgorithm	
	 * 	
	 * @return javax.swing.JPanel	
	 */
	private JPanel getJPanelAlgorithm() {
		if (jPanelAlgorithm == null) {
			jLabelDt = new JLabel();
			jLabelDt.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelDt.setBounds(new Rectangle(60, 211, 122, 25));
			jLabelDt.setHorizontalAlignment(SwingConstants.LEFT);
			jLabelDt.setText("Integration dt, s");
			jLabel1111 = new JLabel();
			jLabel1111.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel1111.setBounds(new Rectangle(60, 111, 106, 25));
			jLabel1111.setText("Model of tumble");
			jLabel111 = new JLabel();
			jLabel111.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel111.setText("Model of run");
			jLabel111.setBounds(new Rectangle(60, 61, 91, 25));
			jLabel115 = new JLabel();
			jLabel115.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel115.setBounds(new Rectangle(60, 161, 122, 25));
			jLabel115.setText("Model of motor");
			jPanelAlgorithm = new JPanel();
			jPanelAlgorithm.setLayout(null);
			jPanelAlgorithm.add(jLabel111, null);
			jPanelAlgorithm.add(getJComboBoxModelRun(), null);
			jPanelAlgorithm.add(jLabel1111, null);
			jPanelAlgorithm.add(getJComboBoxModelTumble(), null);
			jPanelAlgorithm.add(getJComboBoxModelMotor(), null);
			jPanelAlgorithm.add(jLabel115, null);
			jPanelAlgorithm.add(jLabelDt, null);
			jPanelAlgorithm.add(getJTextFieldDt(), null);
		}
		return jPanelAlgorithm;
	}

	/**
	 * This method initializes jTextFieldDt	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldDt() {
		if (jTextFieldDt == null) {
			jTextFieldDt = new JTextField();
			jTextFieldDt.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jTextFieldDt.setBounds(new Rectangle(195, 211, 37, 25));
			jTextFieldDt.setHorizontalAlignment(JTextField.RIGHT);
			jTextFieldDt.setToolTipText("Time step for updating the cell state");
			jTextFieldDt.setText("0.01");
		}
		return jTextFieldDt;
	}

	/**
	 * This method initializes jTextFieldSeed	
	 * 	
	 * @return javax.swing.JTextField	
	 */
	private JTextField getJTextFieldSeed() {
		if (jTextFieldSeed == null) {
			jTextFieldSeed = new JTextField();
			jTextFieldSeed.setFont(new Font("Tahoma", Font.PLAIN, 14));
			jTextFieldSeed.setHorizontalAlignment(SwingConstants.RIGHT);
			jTextFieldSeed.setBounds(new Rectangle(657, 82, 37, 20));
			jTextFieldSeed.setToolTipText("Put a new number to have different simulation");
			jTextFieldSeed.setText("3");
		}
		return jTextFieldSeed;
	}

	/**
	 * This method initializes jComboBoxBoundary	
	 * 	
	 * @return javax.swing.JComboBox	
	 */
	private JComboBox getJComboBoxBoundary() {
		if (jComboBoxBoundary == null) {
			String[] Boundary = {"Reflective", "Periodic"};
			jComboBoxBoundary = new JComboBox(Boundary);
			jComboBoxBoundary.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jComboBoxBoundary.setBounds(new Rectangle(481, 54, 109, 25));
			jComboBoxBoundary.setSelectedIndex(0);
			jComboBoxBoundary.setBackground(Color.white);
			jComboBoxBoundary.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					if(jComboBoxBoundary.getSelectedIndex()==0){
						textAreaConsole.append("\n Walls are solid, cell reflects from wall and makes a tumble.");
					}
					if(jComboBoxBoundary.getSelectedIndex()==1){
						textAreaConsole.append("\n Cell penetrates the wall to appear on the opposite side of arena.");
					}

				}
			});
		}
		return jComboBoxBoundary;
	}
	private ArenaCanvas getCanvas() {
		if (canvas == null) {
			canvas = new ArenaCanvas(Double.parseDouble(jTextFieldIniPosX.getText())/Double.parseDouble(jTextFieldLx.getText()),
									Double.parseDouble(jTextFieldIniPosY.getText())/Double.parseDouble(jTextFieldLy.getText()),
									Double.parseDouble(jTextFieldLx.getText()), Double.parseDouble(jTextFieldLy.getText()),
									Integer.parseInt(jTextFieldNcells.getText()));
			canvas.setLocation(35, 45);
			canvas.setForeground(Color.RED);
			canvas.setBackground(Color.CYAN);
		}
		return canvas;
	}
	private void setCanvas(){
		canvas.setCanvas(Double.parseDouble(jTextFieldIniPosX.getText())/Double.parseDouble(jTextFieldLx.getText()),
							Double.parseDouble(jTextFieldIniPosY.getText())/Double.parseDouble(jTextFieldLy.getText()),
							Double.parseDouble(jTextFieldLx.getText()), Double.parseDouble(jTextFieldLy.getText()),
							Integer.parseInt(jTextFieldNcells.getText()));
	}
	private JPanel getPanelArena() {
		if (panelArena == null) {
			panelArena = new JPanel();
			panelArena.setBounds(18, 255, 723, 62);
			panelArena.setLayout(null);
			jLabel14 = new JLabel();
			jLabel14.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabel14.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel14.setBounds(0, 0, 186, 20);
			panelArena.add(jLabel14);
			jLabel14.setText("Ini position (X,Y), mm");
			jLabel14.setToolTipText("Initial position of all cells. Values in (0,Lx) and (0,Ly), in mm.");
			panelArena.add(getJComboBoxIniPosition());
			panelArena.add(getJTextFieldIniPosX());
			panelArena.add(getJTextFieldIniPosY());
			jLabelIniOri = new JLabel();
			jLabelIniOri.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabelIniOri.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelIniOri.setBounds(342, 0, 153, 20);
			panelArena.add(jLabelIniOri);
			jLabelIniOri.setText("Ini orientation");
			panelArena.add(getJComboBoxIniOri());
		}
		return panelArena;
	}
	private JPanel getJpanelCluster() {
		if (jpanelCluster == null) {
			jpanelCluster = new JPanel();
			jpanelCluster.setBounds(536, 159, 239, 134);
			jpanelCluster.setLayout(null);
			jpanelCluster.add(getJTextFieldNtar());
			jpanelCluster.add(getJTextFieldNtsr());
			jLabel1111112 = new JLabel();
			jLabel1111112.setBounds(0, 0, 199, 20);
			jpanelCluster.add(jLabel1111112);
			jLabel1111112.setFont(new Font("Tahoma", Font.PLAIN, 16));
			jLabel1111112.setText("Cluster composition");
			jLabel112 = new JLabel();
			jLabel112.setBounds(0, 24, 37, 18);
			jpanelCluster.add(jLabel112);
			jLabel112.setFont(new Font("Tahoma", Font.BOLD, 15));
			jLabel112.setHorizontalAlignment(SwingConstants.LEADING);
			jLabel112.setText("Tar");
			jLabel1121 = new JLabel();
			jLabel1121.setBounds(85, 24, 37, 18);
			jpanelCluster.add(jLabel1121);
			jLabel1121.setFont(new Font("Tahoma", Font.BOLD, 15));
			jLabel1121.setText("Tsr");
			jLabel1121.setHorizontalAlignment(SwingConstants.LEADING);
			jpanelCluster.add(getJTextFieldMeth());
			jLabel120 = new JLabel();
			jLabel120.setBounds(0, 81, 117, 20);
			jpanelCluster.add(jLabel120);
			jLabel120.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabel120.setText("Ini. meth. (0..8)");
		}
		return jpanelCluster;
	}
	private JPanel getPanelCheckBoxes() {
		if (panelCheckBoxes == null) {
			panelCheckBoxes = new JPanel();
			panelCheckBoxes.setBounds(102, 72, 568, 203);
			panelCheckBoxes.setLayout(null);
			jLabelOutXY = new JLabel();
			jLabelOutXY.setBounds(29, 34, 127, 25);
			panelCheckBoxes.add(jLabelOutXY);
			jLabelOutXY.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelOutXY.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabelOutXY.setText("(X, Y) position");
			panelCheckBoxes.add(getJCheckBoxOutXY());
			jLabelOutOri = new JLabel();
			jLabelOutOri.setBounds(48, 69, 108, 25);
			panelCheckBoxes.add(jLabelOutOri);
			jLabelOutOri.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelOutOri.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabelOutOri.setText("Cell orientation");
			panelCheckBoxes.add(getJCheckBoxOutOri());
			jLabelOutCheA = new JLabel();
			jLabelOutCheA.setBounds(65, 106, 91, 25);
			panelCheckBoxes.add(jLabelOutCheA);
			jLabelOutCheA.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelOutCheA.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabelOutCheA.setText("CheA-P");
			panelCheckBoxes.add(getJCheckBoxOutCheA());
			jLabelOutCheY = new JLabel();
			jLabelOutCheY.setBounds(65, 143, 91, 25);
			panelCheckBoxes.add(jLabelOutCheY);
			jLabelOutCheY.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelOutCheY.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabelOutCheY.setText("CheY-P");
			panelCheckBoxes.add(getJCheckBoxOutCheY());
			jLabelOutMeth = new JLabel();
			jLabelOutMeth.setBounds(0, 178, 156, 25);
			panelCheckBoxes.add(jLabelOutMeth);
			jLabelOutMeth.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelOutMeth.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabelOutMeth.setText("Receptor methylation");
			panelCheckBoxes.add(getJCheckBoxOutMeth());
			jLabelOutMb = new JLabel();
			jLabelOutMb.setBounds(264, 34, 109, 25);
			panelCheckBoxes.add(jLabelOutMb);
			jLabelOutMb.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelOutMb.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabelOutMb.setText("CCW motor bias");
			panelCheckBoxes.add(getJCheckBoxOutMb());
			jLabelOutRun = new JLabel();
			jLabelOutRun.setBounds(234, 140, 144, 25);
			panelCheckBoxes.add(jLabelOutRun);
			jLabelOutRun.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelOutRun.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabelOutRun.setToolTipText("Value 1 if cell is running, 0 if tumbling at given time t");
			jLabelOutRun.setText("Running or tumbling");
			panelCheckBoxes.add(getJCheckBoxOutRun());
			jLabelOutAttr = new JLabel();
			jLabelOutAttr.setBounds(184, 177, 191, 25);
			panelCheckBoxes.add(jLabelOutAttr);
			jLabelOutAttr.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelOutAttr.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabelOutAttr.setToolTipText("Current attractant concentration experienced by the cell at given time and position");
			jLabelOutAttr.setText("[attr.] sensed by cell, mM");
			panelCheckBoxes.add(getJCheckBoxOutAttr());
			jLabelOutNmot = new JLabel();
			jLabelOutNmot.setBounds(264, 106, 109, 25);
			panelCheckBoxes.add(jLabelOutNmot);
			jLabelOutNmot.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelOutNmot.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabelOutNmot.setToolTipText("Number of motors rotating CW at given time t");
			jLabelOutNmot.setText("N motors CW");
			panelCheckBoxes.add(getJCheckBoxOutNmot());
			jLabelOutSwitchFreq = new JLabel();
			jLabelOutSwitchFreq.setBounds(217, 69, 156, 25);
			panelCheckBoxes.add(jLabelOutSwitchFreq);
			jLabelOutSwitchFreq.setFont(new Font("Tahoma", Font.PLAIN, 15));
			jLabelOutSwitchFreq.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabelOutSwitchFreq.setText("Freq.(CCW->CW), 1/s");
			panelCheckBoxes.add(getJCheckBoxOutSwitchFreq());
			jLabelCheckAll = new JLabel();
			jLabelCheckAll.setBounds(453, 143, 80, 25);
			panelCheckBoxes.add(jLabelCheckAll);
			jLabelCheckAll.setFont(new Font("Tahoma", Font.BOLD, 15));
			jLabelCheckAll.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabelCheckAll.setText("check all");
			panelCheckBoxes.add(getJCheckBoxCheckAll());
			jLabelUncheckAll = new JLabel();
			jLabelUncheckAll.setBounds(430, 178, 103, 25);
			panelCheckBoxes.add(jLabelUncheckAll);
			jLabelUncheckAll.setFont(new Font("Tahoma", Font.BOLD, 15));
			jLabelUncheckAll.setHorizontalAlignment(SwingConstants.RIGHT);
			jLabelUncheckAll.setText("uncheck all");
			panelCheckBoxes.add(getJCheckBoxUncheckAll());
			
			JLabel lblWriteWhat = new JLabel("Write in output files");
			lblWriteWhat.setBounds(48, 0, 226, 25);
			panelCheckBoxes.add(lblWriteWhat);
			lblWriteWhat.setFont(new Font("Tahoma", Font.BOLD, 15));
			
			jCheckBoxOutPon = new JCheckBox("");
			jCheckBoxOutPon.setBounds(543, 34, 25, 25);
			panelCheckBoxes.add(jCheckBoxOutPon);
			jCheckBoxOutPon.setHorizontalAlignment(SwingConstants.TRAILING);
			
			JLabel lblProbPon = new JLabel("Probability P_on");
			lblProbPon.setBounds(430, 36, 103, 20);
			panelCheckBoxes.add(lblProbPon);
			lblProbPon.setToolTipText("Probability that receptor cluster is ON");
			lblProbPon.setHorizontalAlignment(SwingConstants.RIGHT);
			lblProbPon.setFont(new Font("Tahoma", Font.PLAIN, 15));
		}
		return panelCheckBoxes;
	}
}  //  @jve:decl-index=0:visual-constraint="10,16"
