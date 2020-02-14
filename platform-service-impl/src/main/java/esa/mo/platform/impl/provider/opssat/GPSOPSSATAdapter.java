/* ----------------------------------------------------------------------------
 * Copyright (C) 2015      European Space Agency
 *                         European Space Operations Centre
 *                         Darmstadt
 *                         Germany
 * ----------------------------------------------------------------------------
 * System                : ESA NanoSat MO Framework
 * ----------------------------------------------------------------------------
 * Licensed under the European Space Agency Public License, Version 2.0
 * You may not use this file except in compliance with the License.
 *
 * Except as expressly set forth in this License, the Software is provided to
 * You on an "as is" basis and without warranties of any kind, including without
 * limitation merchantability, fitness for a particular purpose, absence of
 * defects or errors, accuracy or non-infringement of intellectual property rights.
 * 
 * See the License for the specific language governing permissions and
 * limitations under the License. 
 * ----------------------------------------------------------------------------
 */
package esa.mo.platform.impl.provider.opssat;

import esa.mo.nanomind.impl.util.NanomindServicesConsumer;
import esa.mo.platform.impl.provider.gen.GPSNMEAonlyAdapter;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.ccsds.moims.mo.mal.MALException;
import esa.opssat.nanomind.opssat_pf.gps.consumer.GPSAdapter;
import org.ccsds.moims.mo.mal.MALInteractionException;
import org.ccsds.moims.mo.mal.structures.Blob;
import org.ccsds.moims.mo.platform.gps.structures.TwoLineElementSet;

/**
 *
 * @author Cesar Coelho
 */
public class GPSOPSSATAdapter extends GPSNMEAonlyAdapter
{

  private final NanomindServicesConsumer obcServicesConsumer;

  public GPSOPSSATAdapter(NanomindServicesConsumer obcServicesConsumer)
  {
    this.obcServicesConsumer = obcServicesConsumer;
  }

  @Override
  public String getNMEASentence(String identifier) throws IOException
  {
    GPSHandler gpsHandler = new GPSHandler();
    try {
      obcServicesConsumer.getGPSNanomindService().getGPSNanomindStub().getGPSData(new Blob(
          identifier.getBytes()), gpsHandler);
    } catch (MALInteractionException | MALException ex) {
      Logger.getLogger(GPSOPSSATAdapter.class.getName()).log(Level.SEVERE, null, ex);
      throw new IOException("Error when retrieving GPS NMEA response from Nanomind", ex);
    }
    return gpsHandler.response;
  }

  @Override
  public boolean isUnitAvailable()
  {
    return true;
  }

  @Override
  public TwoLineElementSet getTLE()
  {
    throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
  }

  private class GPSHandler extends GPSAdapter
  {

    String response = "";

    @Override
    public void getGPSDataResponseReceived(
        org.ccsds.moims.mo.mal.transport.MALMessageHeader msgHeader,
        org.ccsds.moims.mo.mal.structures.Blob data, java.util.Map qosProperties)
    {
      try {
        response = Arrays.toString(data.getValue());
      } catch (MALException ex) {
        Logger.getLogger(GPSHandler.class.getName()).log(Level.SEVERE, null, ex);
      }
    }
  }
}
