package jguis;

import javax.swing.JTable;
import javax.swing.event.TableModelEvent;

public interface TableDialogListener{

	public void tableDataChanged(JTable table,TableModelEvent e,Object[][] tabledata);
	
}
